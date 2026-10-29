#!/usr/bin/env python3
"""One analysis rule per invocation; Snakemake owns dependencies and scheduling."""
from __future__ import annotations
import argparse, csv, gzip, json, math, re, subprocess, sys, sqlite3, statistics
from collections import Counter, defaultdict
from pathlib import Path
import yaml
from design import load_design

RC = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def op(path, mode='rt'):
    return gzip.open(path, mode) if str(path).endswith('.gz') else open(path, mode)

def fastq(path, allow_empty=False):
    with op(path) as f:
        while True:
            h=f.readline()
            if not h: return
            s=f.readline().rstrip(); plus=f.readline(); q=f.readline().rstrip()
            if not h.startswith('@') or not plus.startswith('+') or len(s)!=len(q) or (not s and not allow_empty):
                raise ValueError(f'Malformed FASTQ in {path}: {h[:80]}')
            yield h.rstrip()[1:], s, q

def js(path, data):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n')

def tsv(path, fields, rows):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with op(path,'wt') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

def readtsv(path):
    with op(path) as f:
        yield from csv.DictReader((l for l in f if not l.startswith('#')), delimiter='\t')

def find(root, pattern):
    paths=list(Path(root).rglob(pattern))
    if len(paths)!=1: raise RuntimeError(f'Expected one {pattern} in {root}, found {paths}')
    return paths[0]

def run(cmd):
    print('COMMAND', ' '.join(map(str,cmd)),flush=True)
    subprocess.run(list(map(str,cmd)),check=True)

class Context:
    def __init__(self,config):
        self.config_path=config
        self.c=yaml.safe_load(Path(config).read_text())
        root=Path(config).resolve().parent.parent
        self.sample_rows,self.contrasts=load_design(root/self.c['samples_tsv'],root/self.c['contrasts_tsv'])
        self.samples={r['sample_id']:r['fastq'] for r in self.sample_rows}
        self.groups={r['sample_id']:r.get('group','') for r in self.sample_rows}
        self.batches={r['sample_id']:r.get('batch','') for r in self.sample_rows}
        self.o=Path(self.c['outputs_root']); self.t=Path(self.c['scratch_root']);self.p=self.o/'_project'
    def out(self,s,stage):
        p=self.o/s/stage;p.mkdir(parents=True,exist_ok=True);return p

class QC:
    def __init__(self):
        self.lengths=Counter();self.n=0;self.bases=0;self.qsum=0;self.q20=0;self.q30=0;self.headq=[]
        self.q20map=bytes(int(i>=53) for i in range(256));self.q30map=bytes(int(i>=63) for i in range(256))
    def add(self,h,s,q):
        b=q.encode('ascii');l=len(s); self.n+=1;self.lengths[l]+=1;self.bases+=l
        self.qsum+=sum(b)-33*l;self.q20+=sum(b.translate(self.q20map));self.q30+=sum(b.translate(self.q30map))
        m=re.search(r'QAverage=([\d.]+)',h)
        if m:self.headq.append(float(m.group(1)))
    def data(self):
        cum=0;median=0;n50=0;half=(self.n+1)//2
        for l,n in sorted(self.lengths.items()):
            cum+=n
            if cum>=half:median=l;break
        cum=0
        for l,n in sorted(self.lengths.items(),reverse=True):
            cum+=l*n
            if cum>=self.bases/2:n50=l;break
        return dict(reads=self.n,bases=self.bases,min_length=min(self.lengths,default=0),max_length=max(self.lengths,default=0),median_length=median,n50=n50,mean_length=self.bases/max(1,self.n),mean_base_quality=self.qsum/max(1,self.bases),q20_base_fraction=self.q20/max(1,self.bases),q30_base_fraction=self.q30/max(1,self.bases),header_QAverage_mean=statistics.mean(self.headq) if self.headq else None,quality_encoding='Phred+33; header QAverage reported separately')

def qc(x,s,threads):
    q=QC()
    for h,seq,qual in fastq(x.samples[s]):q.add(h,seq,qual)
    if not q.n:raise RuntimeError('Empty FASTQ')
    d=x.out(s,'01_raw_qc');js(d/'qc.json',q.data());tsv(d/'lengths.tsv',['length','reads'],({'length':l,'reads':n} for l,n in sorted(q.lengths.items())))

def reference(x,s,threads):
    d=x.t/'reference';d.mkdir(parents=True,exist_ok=True)
    g=d/'genes.gtf'
    with op(x.c['reference_gtf']) as src,g.open('w') as dst:
        for line in src:dst.write(line)
    transcripts=defaultdict(list)
    with g.open() as f:
        for line in f:
            if line.startswith('#'):continue
            a=line.rstrip().split('\t')
            if len(a)!=9 or a[2]!='exon':continue
            tid=re.search(r'transcript_id "([^"]+)"',a[8])
            if tid:transcripts[(a[0],a[6],tid[1])].append((int(a[3])-1,int(a[4])))
    with (d/'junctions.bed').open('w') as f:
        for (chrom,strand,tid),exons in sorted(transcripts.items()):
            exons=sorted(set(exons));start=exons[0][0];end=exons[-1][1]
            f.write('\t'.join(map(str,[chrom,start,end,tid,0,strand,start,end,0,len(exons),','.join(str(b-a) for a,b in exons)+',',','.join(str(a-start) for a,b in exons)+',']))+'\n')
    run(['minimap2','-x','splice','-t',str(x.c.get('threads',8)),'-d',d/'genome.splice.mmi',x.c['reference_fasta']])

def orient(x,s,threads):
    d=x.out(s,'02_trim');w=x.t/s;w.mkdir(parents=True,exist_ok=True)
    stats=Counter();quality=QC();rows_path=w/'trim_metadata.tsv.gz';poly=re.compile('T{'+str(x.c['poly_t_min'])+',}')
    with gzip.open(w/'insert.fastq.gz','wt',compresslevel=1) as f,gzip.open(rows_path,'wt',compresslevel=1) as meta:
        writer=csv.writer(meta,delimiter='\t',lineterminator='\n');writer.writerow(['read_id','original_read_id','sample_id','polyA_supported','polyA_length','orientation'])
        for h,seq,q in fastq(w/'adapters.fastq.gz', allow_empty=True):
            stats['dual_adapter_reads']+=1
            m=poly.search(seq[:int(x.c['poly_t_window'])])
            supported=m is not None
            if supported: seq=seq[m.start():];q=q[m.start():]
            # Cutadapt yields primer-left, poly(T)-left cDNA: reverse to mRNA sense.
            seq=seq.translate(RC)[::-1];q=q[::-1]
            if len(seq)<int(x.c['min_insert_length']):stats['short_insert_reads']+=1;continue
            rid=s+'|'+h.split()[0]
            f.write('@'+rid+'\n'+seq+'\n+\n'+q+'\n')
            writer.writerow([rid,h.split()[0],s,int(supported),len(m[0]) if supported else 0,'reverse_input' if h.endswith(' rc') else 'forward_input'])
            stats['retained_reads']+=1;stats['polyA_supported_reads' if supported else 'polyA_uncertain_reads']+=1;quality.add(h,seq,q)
    ca=json.loads((d/'cutadapt.json').read_text());stats['input_reads']=ca['read_counts']['input'];stats['adapter_unassigned_reads']=stats['input_reads']-stats['dual_adapter_reads']
    assert stats['input_reads']==stats['retained_reads']+stats['short_insert_reads']+stats['adapter_unassigned_reads']
    if not stats['retained_reads']:raise RuntimeError('No retained reads after trimming')
    js(d/'trim_summary.json',dict(stats));js(d/'insert_qc.json',quality.data())

def align(x,s,threads):
    import pysam
    d=x.out(s,'03_alignment');w=x.t/s; w.mkdir(parents=True,exist_ok=True)
    bam=d/'aligned.sorted.bam';filtered=d/'primary.q20.bam'
    mapping=max(1,threads-2);sorting=max(1,threads-mapping-1)
    cmd=['minimap2','-ax','splice','-uf','--secondary=no','--MD','--cs=long','--junc-bed',str(x.t/'reference/junctions.bed'),'-t',str(mapping),str(x.t/'reference/genome.splice.mmi'),str(w/'insert.fastq.gz')]
    print('COMMAND',' '.join(cmd),flush=True)
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE)
    sorter=subprocess.Popen(['samtools','sort','-@',str(sorting),'-m','1G','-T',str(w/'sort'),'-o',str(bam),'-'],stdin=p.stdout)
    p.stdout.close(); sr=sorter.wait();mr=p.wait()
    if sr or mr:raise RuntimeError(f'Alignment/sort failed: {mr}/{sr}')
    run(['samtools','index','-@',threads,bam]);run(['samtools','view','-b','-q',x.c['mapq'],'-F','2308','-o',filtered,bam]);run(['samtools','index','-@',threads,filtered])
    c=Counter()
    with pysam.AlignmentFile(bam,'rb') as f:
        for r in f:
            if r.is_secondary or r.is_supplementary:c['secondary_or_supplementary_records']+=1;continue
            c['input_reads']+=1
            if r.is_unmapped:c['unmapped_reads']+=1
            elif r.mapping_quality<int(x.c['mapq']):c['low_mapq_reads']+=1
            else:c['primary_q20_reads']+=1
    assert c['input_reads']==sum(c[k] for k in ['unmapped_reads','low_mapq_reads','primary_q20_reads'])
    js(d/'alignment.json',dict(c))

def sqanti(x,s,threads):
    gtf=find(x.p/'04_catalog/isoquant','*.transcript_models.gtf')
    out=x.p/'05_structure_qc/sqanti3';out.mkdir(parents=True,exist_ok=True)
    run(['sqanti3_qc.py','--isoforms',gtf,'--refGTF',x.t/'reference/genes.gtf','--refFasta',x.t/'reference/genome.fa','--output','catalog','--dir',out,'--cpus',threads,'--report','skip'])

def matrix(path,samples):
    # IsoQuant matrix headers are commented; do not discard them as metadata.
    with op(path) as f:
        lines=[l.rstrip('\n') for l in f if l.strip()]
    header=None;rows={}
    for line in lines:
        a=line.lstrip('#').strip().split('\t')
        if header is None:
            if all(s in a for s in samples):header=a
            continue
        if line.startswith('#') or a[0].startswith('__'):continue
        if len(a)!=len(header):raise RuntimeError(f'Malformed matrix {path}: {line[:100]}')
        values={s:float(a[header.index(s)]) for s in samples}
        if not all(math.isfinite(v) and v>=0 for v in values.values()):raise RuntimeError(f'Invalid abundance in {path}')
        rows[a[0]]=values
    if header is None:raise RuntimeError(f'Missing sample columns in {path}')
    return rows

def model_info(path):
    models={}
    with op(path) as f:
        for l in f:
            if l.startswith('#'):continue
            a=l.rstrip().split('\t')
            if len(a)!=9:continue
            attr=dict(re.findall(r'(\w+) "([^"]*)"',a[8]));tid=attr.get('transcript_id')
            if not tid:continue
            m=models.setdefault(tid,dict(chrom=a[0],strand=a[6],gene_id=attr.get('gene_id','.'),exons=[]))
            if a[2]=='exon':m['exons'].append((int(a[3])-1,int(a[4])))
    for m in models.values():m['exons']=sorted(set(m['exons']))
    return models

def quantify(x,s,threads):
    iq=x.p/'04_catalog/isoquant';qdir=x.p/'06_quantification';qdir.mkdir(parents=True,exist_ok=True)
    counts=matrix(find(iq,'*.discovered_transcript_grouped_file_name_counts.tsv'),x.samples)
    tpm=matrix(find(iq,'*.discovered_transcript_grouped_file_name_tpm.tsv'),x.samples)
    # Grouped discovery counts must reconstruct the native pooled counts.
    pooled={}
    with op(find(iq,'*.discovered_transcript_counts.tsv')) as f:
        for line in f:
            if line.startswith('#') or line.startswith('__'):continue
            a=line.strip().split('\t')
            if len(a)>=2:
                try:pooled[a[0]]=float(a[1])
                except ValueError:continue
    for tid,values in counts.items():
        if tid not in pooled or not math.isclose(sum(values.values()),pooled[tid],abs_tol=1e-6):raise RuntimeError(f'Grouped count mismatch for {tid}')
    gtf=find(iq,'*.transcript_models.gtf');models=model_info(gtf)
    cls=list(readtsv(find(x.p/'05_structure_qc/sqanti3','*classification.txt')))
    assert {r['isoform'] for r in cls}==set(models), 'SQANTI/IsoQuant model IDs differ'
    catalog=[];keep=set();cc=Counter();genes=Counter()
    for r in cls:
        tid=r['isoform'];cat=r['structural_category'];m=models[tid]
        minimum=1 if cat=='full-splice_match' else 2 if cat=='novel_in_catalog' else 3
        if len(m['exons'])==1 and cat!='full-splice_match':minimum=max(minimum,5)
        c=counts.get(tid,{sid:0. for sid in x.samples});reasons=[]
        if max(c.values())<minimum:reasons.append(f'support_lt_{minimum}')
        if r.get('RTS_stage','').lower() in ['true','1']:reasons.append('RT_switching')
        if r.get('all_canonical','').lower() in ['non_canonical','false','0']:reasons.append('noncanonical_splice')
        try:
            if float(r.get('perc_A_downstream_TTS','0'))>=60:reasons.append('intrapriming')
        except ValueError:pass
        passed=not reasons
        if passed:keep.add(tid)
        cc[cat]+=1;genes[m['gene_id']]+=1
        catalog.append(dict(transcript_id=tid,gene_id=m['gene_id'],gene_name=r.get('associated_gene','.'),associated_transcript=r.get('associated_transcript','.'),structural_category=cat,chrom=m['chrom'],strand=m['strand'],length=r.get('length',''),exons=len(m['exons']),filter_status='PASS' if passed else 'FAIL',filter_reason=';'.join(reasons) or 'PASS',**c))
    if not catalog:raise RuntimeError('No models classified')
    tsv(qdir/'catalog.tsv',list(catalog[0]),catalog)
    with op(gtf) as src,(qdir/'filtered.gtf').open('w') as dst:
        for l in src:
            match=re.search(r'transcript_id "([^"]+)"',l)
            if match and match[1] in keep:dst.write(l)
    for name,mat in [('isoform_counts',counts),('isoform_tpm',tpm)]:
        tsv(qdir/(name+'.tsv'),['transcript_id',*x.samples],({'transcript_id':tid,**mat.get(tid,{s:0. for s in x.samples})} for tid in sorted(keep)))
    # Retain exact native reference matrices, independent of model filtering.
    refsummary={}
    for feature in ['gene','transcript']:
        for unit in ['counts','tpm']:
            native=matrix(find(iq,f'*.{feature}_grouped_file_name_{unit}.tsv'),x.samples)
            tsv(qdir/f'reference_{feature}_{unit}.tsv',[feature+'_id',*x.samples],({feature+'_id':tid,**values} for tid,values in sorted(native.items())))
            if unit=='counts':refsummary[feature]={sid:sum(v[sid]>0 for v in native.values()) for sid in x.samples}
    totals={sid:sum(v[sid] for v in counts.values()) for sid in x.samples}
    retained={sid:sum(counts.get(tid,{}).get(sid,0) for tid in keep) for sid in x.samples}
    usages=[]
    for sid in x.samples:
        gt=Counter()
        for tid in keep:gt[models[tid]['gene_id']]+=counts.get(tid,{}).get(sid,0)
        for tid in sorted(keep):
            gene=models[tid]['gene_id'];n=counts.get(tid,{}).get(sid,0)
            usages.append(dict(sample_id=sid,transcript_id=tid,gene_id=gene,reads=n,gene_reads=gt[gene],isoform_usage=n/gt[gene] if gt[gene] else 'NA'))
    tsv(qdir/'isoform_usage.tsv',['sample_id','transcript_id','gene_id','reads','gene_reads','isoform_usage'],usages)
    for sid in x.samples:
        tsv(x.out(sid,'06_quantification')/'isoform_counts.tsv',['transcript_id','reads'],({'transcript_id':tid,'reads':counts.get(tid,{}).get(sid,0)} for tid in sorted(keep)))
    js(qdir/'summary.json',dict(catalog_models=len(catalog),retained_models=len(keep),structural_categories=dict(cc),model_assigned_reads=totals,retained_model_reads=retained,reference_detected=refsummary,counting_unit='unique-consistent reads; not UMI molecules',reference='GRCm39 / supplied filtered GENCODE M33',samples=list(x.samples)))

def end_clusters(points,window):
    groups=[]
    for p in sorted(points,key=lambda r:r[0]):
        if not groups or p[0]-groups[-1][0][0]>window:groups.append([])
        groups[-1].append(p)
    return groups

def events(x,s,threads):
    import pysam
    d=x.p/'07_events';d.mkdir(parents=True,exist_ok=True);q=x.p/'06_quantification'
    counts=matrix(q/'isoform_counts.tsv',x.samples);tpm=matrix(q/'isoform_tpm.tsv',x.samples);models=model_info(q/'filtered.gtf')
    events_rows=[]
    if models:
        expression=d/'suppa_expression.tsv'
        with expression.open('w') as f:
            f.write('\t'.join(x.samples)+'\n')
            for tid in sorted(models):f.write(tid+'\t'+'\t'.join(str(tpm.get(tid,{}).get(sid,0.)) for sid in x.samples)+'\n')
        run(['suppa.py','generateEvents','-i',q/'filtered.gtf','-o',d/'events','-f','ioe','-e','SE','SS','MX','RI','FL'])
        for ioe in sorted(d.glob('events*.ioe')):
            if sum(1 for _ in ioe.open())<=1:continue
            run(['suppa.py','psiPerEvent','-i',ioe,'-e',expression,'-o',ioe.with_suffix('')])
            psipath=ioe.with_suffix('.psi');psimat={}
            with psipath.open() as f:
                header=f.readline().strip().split('\t')
                for line in f:
                    a=line.strip().split('\t');psimat[a[0]]=dict(zip(header,a[1:]))
            for r in readtsv(ioe):
                eid=r['event_id'];denids=r['total_transcripts'].split(',')
                for sid in x.samples:
                    support=sum(counts.get(tid,{}).get(sid,0) for tid in denids)
                    value=psimat.get(eid,{}).get(sid,'nan');psi=float(value) if value.lower() not in ['na','nan'] else None
                    events_rows.append(dict(event_id=eid,sample_id=sid,gene_id=r['gene_id'],event_type=eid.split(';')[-1].split(':')[0],support_reads=support,psi=psi if support>=x.c['event_min_reads'] and psi is not None else 'NA'))
    tsv(d/'psi.tsv',['event_id','sample_id','gene_id','event_type','support_reads','psi'],events_rows)
    # Disk-backed read/model/tail join keeps full data memory bounded.
    db=x.t/'events.sqlite';db.parent.mkdir(parents=True,exist_ok=True)
    if db.exists():db.unlink()
    conn=sqlite3.connect(db);conn.execute('PRAGMA journal_mode=OFF');conn.execute('PRAGMA synchronous=OFF')
    conn.execute('CREATE TABLE reads (rid TEXT PRIMARY KEY, tid TEXT, multiple INTEGER DEFAULT 0, tail INTEGER DEFAULT 0)')
    r2t=find(x.p/'04_catalog/isoquant','*.read_assignments.tsv.gz')
    batch=[]
    for r in readtsv(r2t):
        tid=r['isoform_id']
        if tid in models:
            batch.append((r['read_id'],tid))
            if len(batch)>=50000:
                conn.executemany('INSERT INTO reads(rid,tid) VALUES(?,?) ON CONFLICT(rid) DO UPDATE SET multiple=1',batch);batch=[]
    conn.executemany('INSERT INTO reads(rid,tid) VALUES(?,?) ON CONFLICT(rid) DO UPDATE SET multiple=1',batch);conn.commit()
    tail_updates=[]
    for sid in x.samples:
        for r in readtsv(x.t/sid/'trim_metadata.tsv.gz'):
            if r['polyA_supported']=='1':tail_updates.append((r['read_id'],))
            if len(tail_updates)>=50000:conn.executemany('UPDATE reads SET tail=1 WHERE rid=?',tail_updates);tail_updates=[]
    conn.executemany('UPDATE reads SET tail=1 WHERE rid=?',tail_updates);conn.commit()
    terminal=[];eligible={}
    for sid in x.samples:
        points=defaultdict(list);n=0
        with pysam.AlignmentFile(x.o/sid/'03_alignment/primary.q20.bam','rb') as bam:
            for read in bam:
                row=conn.execute('SELECT tid,tail,multiple FROM reads WHERE rid=?',(read.query_name,)).fetchone()
                if not row or row[2]:continue
                tid,tail,_=row;m=models[tid]
                if ('-' if read.is_reverse else '+')!=m['strand']:continue
                if read.query_alignment_length/max(1,read.query_length)<0.8:continue
                n+=1
                start=read.reference_start;end=read.reference_end-1
                for kind,pos in [('5prime',end if m['strand']=='-' else start),('3prime',start if m['strand']=='-' else end)]:
                    if kind=='3prime' and not tail:continue
                    points[(m['gene_id'],m['chrom'],m['strand'],kind)].append((pos,tid))
        eligible[sid]=n
        for (gene,chrom,strand,kind),pts in sorted(points.items()):
            for group in end_clusters(pts,int(x.c['end_window'])):
                if len(group)<int(x.c['end_min_reads']):continue
                terminal.append(dict(sample_id=sid,gene_id=gene,chrom=chrom,strand=strand,end_type=kind,position_1based=int(statistics.median(p[0] for p in group))+1,support_reads=len(group),transcripts=','.join(sorted({p[1] for p in group})),evidence='polyA-supported read ends' if kind=='3prime' else 'aligned read starts; TSS unvalidated'))
    conn.close()
    tsv(d/'terminal_candidates.tsv',['sample_id','gene_id','chrom','strand','end_type','position_1based','support_reads','transcripts','evidence'],terminal)
    js(d/'events.json',dict(psi_rows=len(events_rows),evaluable_psi_rows=sum(r['psi']!='NA' for r in events_rows),terminal_candidates=len(terminal),eligible_terminal_reads=eligible,status='complete',interpretation='descriptive only; no biological p-values; aligned starts are not validated TSS; candidate ends are not validated APA sites'))

def complete(x,s,threads):
    report=x.p/'09_report'
    for sid in x.samples:
        raw=json.loads((x.o/sid/'01_raw_qc/qc.json').read_text());trim=json.loads((x.o/sid/'02_trim/trim_summary.json').read_text());al=json.loads((x.o/sid/'03_alignment/alignment.json').read_text())
        if raw['reads']!=trim['input_reads'] or trim['retained_reads']!=al['input_reads']:raise RuntimeError(f'Read-flow mismatch for {sid}')
    run([sys.executable,'pipeline/scripts/validate_delivery_report.py','--report',report/'report.html','--manifest',report/'delivery_manifest.tsv','--require-section','1 项目概览','--require-expandable-figures','--require-print-css','--json'])
    # Completion reflects validated required artifacts, never merely a successful PID.
    for row in readtsv(report/'delivery_manifest.tsv'):
        p=Path(row['target'])
        if not p.is_file() or not p.stat().st_size:raise RuntimeError(f'Missing deliverable: {p}')
    tsv(x.p/'RUN_COMPLETE.tsv',['status','samples','report'],[dict(status='completed',samples=len(x.samples),report=str(report/'report.html'))])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['qc','reference','orient','align','sqanti','quantify','events','complete']);p.add_argument('--config',required=True);p.add_argument('--sample');p.add_argument('--threads',type=int,default=8)
    a=p.parse_args();globals()[a.stage](Context(a.config),a.sample,a.threads)

if __name__=='__main__':main()
