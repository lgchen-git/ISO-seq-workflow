#!/usr/bin/env python3
"""Generate the sole offline customer report from completed analytical tables."""
from __future__ import annotations
import argparse, base64, html, json, subprocess, datetime, importlib.metadata, os, hashlib, colorsys
from pathlib import Path
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from analyze import Context, readtsv, js, tsv, model_info

COLORS=['#E64B35','#4DBBD5','#00A087','#3C5488','#F39B7F','#8491B4','#91D1C2','#DC0000','#7E6148','#B09C85','#6F99AD','#8A7AA9']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','axes.titleweight':'normal'})

def esc(x):return html.escape(str(x))
def load(p):return json.loads(Path(p).read_text())
def sample_color(sample):
    hue=int.from_bytes(hashlib.sha256(sample.encode('utf-8')).digest()[:4],'big')/2**32
    return colorsys.hsv_to_rgb(hue,.62,.78)
def pvalue_key(row,key):
    try:return float(row.get(key,'1'))
    except (TypeError,ValueError):return 1.0

def table(rows,columns,id,paginate=False):
    attrs=f'id="{id}"'+(' class="paginated-table" data-page-size="10"' if paginate else '')
    text=f'<div class="table-wrap"><table {attrs}><thead><tr>'+''.join('<th>'+esc(label)+'</th>' for key,label in columns)+'</tr></thead><tbody>'
    for r in rows:
        text+='<tr>'+''.join('<td>'+esc(r.get(key,''))+'</td>' for key,label in columns)+'</tr>'
    text+='</tbody></table></div>'
    if paginate:
        text+=f'<div class="pager" data-pager-for="{id}"><input type="search" aria-label="筛选结果" placeholder="筛选当前表"><button data-page-prev>上一页</button><span data-page-status></span><button data-page-next>下一页</button></div>'
    return text

def figure(d,name,fig,title,caption):
    fig.tight_layout();fig.savefig(d/(name+'.png'),dpi=350);fig.savefig(d/(name+'.pdf'));plt.close(fig)
    data=base64.b64encode((d/(name+'.png')).read_bytes()).decode()
    return f'<article class="figure-card expandable" data-expandable="true"><h3>{esc(title)}</h3><button class="figure-expand" aria-label="放大 {esc(title)}"><img src="data:image/png;base64,{data}" alt="{esc(title)}"></button><p class="caption">{esc(caption)}</p></article>'

CSS='''*{box-sizing:border-box}body{margin:0;color:#29313b;background:#f5f6f8;font:15px/1.7 Arial,"Noto Sans SC",sans-serif}.sidebar{position:fixed;left:0;top:0;bottom:0;width:226px;padding:25px 18px;background:#26354b;color:white;overflow:auto}.sidebar a{display:block;padding:10px 0;color:#f1f4f8;text-decoration:none}.sidebar strong{display:block;font-size:19px;margin-bottom:22px}main{margin-left:226px;padding:28px;max-width:1660px}.report-section{background:white;margin-bottom:22px;padding:26px;border:1px solid #e2e6ea;border-radius:8px;scroll-margin-top:18px}h1{font-size:30px;margin:0 0 8px}h2{font-size:23px;margin-top:0}h3{font-size:17px}.subtitle,.caption{color:#627080}.kpis{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:20px 0}.kpi{background:#f6f8fa;border-top:3px solid #3C5488;padding:14px}.kpi span{display:block;font-size:12px;color:#657386}.kpi strong{font-size:23px}.notice{background:#f4f7fa;border-left:3px solid #4DBBD5;padding:12px 16px;margin:14px 0}.grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:20px}.figure-card{border:1px solid #e1e5e9;border-radius:6px;padding:15px;margin:16px 0}.figure-expand{background:none;border:none;padding:0;width:100%;cursor:zoom-in}.figure-expand img{display:block;width:100%;height:auto}.caption{font-size:13px}.table-wrap{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px}th,td{padding:8px 10px;text-align:left;border-bottom:1px solid #e5e8ed}th{background:#f1f4f7;white-space:nowrap}td{max-width:350px;overflow-wrap:anywhere}.pager{display:flex;align-items:center;gap:10px;margin:14px 0}.pager input{flex:1;min-width:80px}input,button{font:inherit}input,.pager button,#close-modal{padding:7px 10px;border:1px solid #cbd3dd;border-radius:4px;background:white}details{margin:16px 0;padding:12px;border:1px solid #dde3e9}summary{cursor:pointer;font-weight:bold}.footer{margin-left:226px;padding:20px 28px;background:#26354b;color:white}#modal{position:fixed;inset:0;z-index:10;background:#000b;padding:45px 20px;overflow:auto}#modal[hidden]{display:none}#modal img{display:block;max-width:95%;margin:auto;background:white}#close-modal{position:fixed;right:25px;top:10px}button:focus-visible,a:focus-visible,input:focus-visible{outline:3px solid #4DBBD5}@media(max-width:1000px){.kpis{grid-template-columns:repeat(2,1fr)}.grid{grid-template-columns:1fr}}@media(max-width:740px){.sidebar{position:static;width:auto}main,.footer{margin-left:0}main{padding:12px}.report-section{padding:17px}h1{font-size:24px}.pager{flex-wrap:wrap}}@media print{body{background:white}.sidebar,#modal,.pager{display:none!important}main,.footer{margin:0;padding:0}.report-section{border:0;padding:12px 0;box-shadow:none}tr[hidden]{display:table-row!important}.figure-card{break-inside:avoid}.grid{display:block}details{display:block}details>*{display:block!important}.footer{color:#333;background:white}.kpis{grid-template-columns:repeat(4,1fr)}}'''
JS='''const modal=document.getElementById('modal');function closeModal(){modal.hidden=true;modal.querySelector('img').removeAttribute('src')}document.querySelectorAll('.figure-expand').forEach(b=>b.addEventListener('click',()=>{modal.querySelector('img').src=b.querySelector('img').src;modal.querySelector('img').alt=b.querySelector('img').alt;modal.hidden=false;document.getElementById('close-modal').focus()}));document.getElementById('close-modal').addEventListener('click',closeModal);modal.addEventListener('click',e=>{if(e.target===modal)closeModal()});document.addEventListener('keydown',e=>{if(e.key==='Escape')closeModal()});document.querySelectorAll('.paginated-table').forEach(t=>{const rows=[...t.tBodies[0].rows],p=document.querySelector('[data-pager-for="'+t.id+'"]'),size=10;let page=0;function render(){const text=p.querySelector('input').value.toLowerCase(),filtered=rows.filter(r=>r.textContent.toLowerCase().includes(text)),pages=Math.max(1,Math.ceil(filtered.length/size));page=Math.max(0,Math.min(page,pages-1));rows.forEach(r=>r.hidden=true);filtered.slice(page*size,(page+1)*size).forEach(r=>r.hidden=false);p.querySelector('[data-page-status]').textContent=(page+1)+' / '+pages+' 页 · '+filtered.length+' 条';p.querySelector('[data-page-prev]').disabled=page===0;p.querySelector('[data-page-next]').disabled=page>=pages-1}p.querySelector('input').addEventListener('input',()=>{page=0;render()});p.querySelector('[data-page-prev]').addEventListener('click',()=>{page--;render()});p.querySelector('[data-page-next]').addEventListener('click',()=>{page++;render()});render()});let detailStates=[];window.addEventListener('beforeprint',()=>{detailStates=[...document.querySelectorAll('details')].map(d=>[d,d.open]);detailStates.forEach(([d])=>d.open=true)});window.addEventListener('afterprint',()=>detailStates.forEach(([d,s])=>d.open=s));'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',required=True);a=ap.parse_args();x=Context(a.config)
    d=x.p/'09_report';d.mkdir(parents=True,exist_ok=True);fd=d/'figures';fd.mkdir(exist_ok=True)
    q=x.p/'06_quantification';ed=x.p/'07_events';dd=x.p/'08_differential';summary=load(q/'summary.json');ev=load(ed/'events.json');diff=load(dd/'summary.json');samples=list(x.samples)
    qc={s:load(x.o/s/'01_raw_qc/qc.json') for s in samples};trim={s:load(x.o/s/'02_trim/trim_summary.json') for s in samples};al={s:load(x.o/s/'03_alignment/alignment.json') for s in samples}
    cat=list(readtsv(q/'catalog.tsv'));passing=[r for r in cat if r['filter_status']=='PASS'];psi=list(readtsv(ed/'psi.tsv'));ends=list(readtsv(ed/'terminal_candidates.tsv'))
    sections=[]
    total=sum(v['reads'] for v in qc.values());retained=sum(v.get('retained_reads',0) for v in trim.values());mapped=sum(v.get('primary_q20_reads',0) for v in al.values())
    kpis=[('原始 reads',f'{total:,}'),('进入主分析 reads',f'{mapped:,}'),('共享目录异构体',f"{summary['catalog_models']:,}"),('结构过滤后异构体',f"{summary['retained_models']:,}")]
    intro='<h1>'+esc(x.c['report_title'])+'</h1><p class="subtitle">'+esc(x.c['project_id'])+' · Mus musculus · GRCm39 / GENCODE M33 · bulk 长读长 cDNA</p><div class="kpis">'+''.join(f'<div class="kpi"><span>{esc(k)}</span><strong>{esc(v)}</strong></div>' for k,v in kpis)+'</div>'
    intro+=f'<p>{len(samples)} 个独立生物学样本共享转录本目录，分别保留表达与读段支持。双端接头处理后保留 {retained:,} 条 reads；MAPQ≥20 的主比对用于定量和结构分析。各阶段数量及分母见下方质控表。</p>'
    if diff['stats_enabled']:
        intro+='<div class="notice">本报告以 reads 为计数单位，无 UMI 去重。差异分析以样本表中的独立生物学重复为单位；PCR reads 不等同于独立 RNA 分子。</div>'
    else:
        intro+='<div class="notice">本报告以 reads 为计数单位，无 UMI 去重。未提供有效比较表，因此只输出描述性结果，不计算差异表达或 DTU 的 p 值/FDR；PCR reads 不等同于独立 RNA 分子。</div>'
    sections.append(('overview','1 项目概览',intro))
    qcrows=[]
    for s in samples:
        v=qc[s];t=trim[s];b=al[s]
        qcrows.append(dict(sample=s,group=x.groups[s] or '未指定',batch=x.batches[s] or '未指定',raw=f"{v['reads']:,}",median=v['median_length'],n50=v['n50'],q30=f"{v['q30_base_fraction']:.2%}",header=f"{v['header_QAverage_mean']:.2f}" if v['header_QAverage_mean'] is not None else 'NA',retained=f"{t.get('retained_reads',0):,}",tail=f"{t.get('polyA_supported_reads',0)/max(1,t.get('retained_reads',0)):.2%}",mapped=f"{b.get('primary_q20_reads',0):,}",maprate=f"{b.get('primary_q20_reads',0)/max(1,b.get('input_reads',0)):.2%}"))
    content=table(qcrows,[('sample','样本'),('group','组别'),('batch','Batch'),('raw','原始 reads'),('median','中位读长 bp'),('n50','N50 bp'),('q30','FASTQ Q30 碱基'),('header','头部 QAverage'),('retained','保留 reads'),('tail','保留 reads 尾支持'),('mapped','主比对 Q20 reads'),('maprate','主比对/保留 reads')],'qc-summary',True)
    content+='<p class="caption">Q30 由 FASTQ 的 Phred+33 字符计算，头部 QAverage 为测序系统单独提供的均值，二者分开列示。比对保留率以接头处理后 reads 为分母。</p>'
    fig,ax=plt.subplots(figsize=(7.4,4.4))
    for i,s in enumerate(samples):
        hist=list(readtsv(x.o/s/'01_raw_qc/lengths.tsv'));lengths=np.array([int(r['length']) for r in hist]);n=np.array([int(r['reads']) for r in hist]);bins=np.arange(0,10051,100);values,_=np.histogram(lengths,bins,weights=n)
        ax.plot((bins[:-1]+bins[1:])/2,values/max(1,qc[s]['reads']),label=s,color=sample_color(s))
    ax.set(xlabel='Read length (bp)',ylabel='Fraction of raw reads / 100 bp bin',xlim=(0,10000));ax.legend(frameon=False)
    content+='<div class="grid">'+figure(fd,'raw_lengths',fig,'原始读长分布','每 100 bp 分箱，分母为每个样本全部原始 reads；显示 0–10 kb，超出显示范围的 reads 仍纳入质控统计。')
    fig,ax=plt.subplots(figsize=(7.4,4.4));positions=np.arange(3)
    for i,s in enumerate(samples):
        ax.plot(positions,[1,trim[s].get('retained_reads',0)/qc[s]['reads'],al[s].get('primary_q20_reads',0)/qc[s]['reads']],marker='o',color=sample_color(s),label=s)
    ax.set_xticks(positions,['Raw','Adapter-retained','Primary MAPQ>=20']);ax.set(ylabel='Fraction of raw reads',ylim=(0,1.05));ax.legend(frameon=False)
    content+=figure(fd,'read_retention',fig,'读段处理链','各阶段均除以本样本原始 reads；进入主分析的 reads 不包含未比对、低 MAPQ、次级或补充比对。')+'</div>'
    flow_dir=x.p/'09_report/read_flow'
    flow_fragment=(flow_dir/'reports/read_flow_sankey_section.html').read_text(encoding='utf-8')
    flow_manifest=flow_dir/'reports/sankey_manifest.tsv'
    flow_targets=[]
    for item in readtsv(flow_manifest):
        target=(flow_manifest.parent/item['target']).resolve()
        if target.is_file() and target.stat().st_size:flow_targets.append(target)
    content+='<h3>逐样本守恒读流</h3>'+flow_fragment
    content+='<p class="caption">每个下拉视图是一份独立样本文库，单位为 reads；不同样本视图的根节点不可相加。分支计数来自原始 QC、接头处理和比对阶段的精确摘要，并逐样本验证守恒。</p>'
    sections.append(('qc','2 数据质量与处理链',content))
    preferred=['full-splice_match','incomplete-splice_match','novel_in_catalog','novel_not_in_catalog','genic','genic_intron','antisense','intergenic','fusion'];classes=[c for c in preferred if c in summary['structural_categories']]+sorted(set(summary['structural_categories'])-set(preferred))
    fig,ax=plt.subplots(figsize=(8,4.5));values=[summary['structural_categories'][c] for c in classes];ax.barh(classes[::-1],values[::-1],color=[COLORS[i%len(COLORS)] for i in range(len(classes))][::-1]);ax.set(xlabel='Number of shared transcript models')
    content=figure(fd,'structure_classes',fig,'共享异构体结构分类','计数为 SQANTI3 对全部发现模型的结构分类，包含未通过过滤的模型；“新”指相对本次提供的过滤参考注释，不代表已验证的新功能。')
    top=sorted(passing,key=lambda r:sum(float(r[s]) for s in samples),reverse=True)[:200]
    content+='<h3>支持 reads 最高的 200 个过滤后异构体</h3>'+table(top,[('transcript_id','转录本'),('gene_id','基因'),('structural_category','结构类别'),('exons','外显子'),*[(s,s+' reads') for s in samples]],'catalog-top',True)
    content+='<details><summary>结构过滤规则</summary><p>FSM≥1、NIC≥2、其他≥3 reads；新单外显子至少 5 reads。至少一个样本达到阈值。RT-switching、非经典剪接及下游 A 比例≥60% 的内引发候选被标记并过滤。未过滤的完整目录及过滤原因保留在 catalog.tsv。reads 支持不等于生物学重复或独立分子支持。</p></details>'
    sections.append(('catalog','3 转录本目录与结构质控',content))
    content=table([dict(sample=s,genes=summary['reference_detected']['gene'][s],transcripts=summary['reference_detected']['transcript'][s],modelreads=f"{summary['model_assigned_reads'][s]:,.0f}",filtered=f"{summary['retained_model_reads'][s]:,.0f}") for s in samples],[('sample','样本'),('genes','参考基因 ≥1 read'),('transcripts','参考转录本 ≥1 read'),('modelreads','发现模型分配 reads'),('filtered','过滤后模型 reads')],'quant-qc',True)
    group_names=sorted({group for group in x.groups.values() if group})
    if group_names:
        group_rows=[dict(group=group,sample_count=sum(x.groups[s]==group for s in samples),samples=', '.join(s for s in samples if x.groups[s]==group),batches=', '.join(sorted({x.batches[s] for s in samples if x.groups[s]==group and x.batches[s]})) or '未指定') for group in group_names]
        content+='<h3>组别与生物学重复</h3>'+table(group_rows,[('group','组别'),('sample_count','样本数'),('samples','独立生物学重复'),('batches','Batch')],'group-summary',True)
    content+='<p class="notice">参考注释表达与发现模型表达采用不同分配目录，计数不要求相等，也不混用其分母。基因内使用比例以该样本过滤后模型 reads 为分母。</p>'
    sections.append(('expression','4 表达定量与样本对照',content))
    comparisons=[dict(contrast_id=r['contrast_id'],test_group=r['test_group'],reference_group=r['reference_group'],samples=', '.join(s for s in samples if x.groups[s] in {r['test_group'],r['reference_group']})) for r in x.contrasts]
    de_gene=list(readtsv(dd/'gene_de.tsv'));de_tx=list(readtsv(dd/'transcript_de.tsv'));dtu_gene=list(readtsv(dd/'dtu_genes.tsv'));dtu_tx=list(readtsv(dd/'dtu_transcripts.tsv'))
    content=''
    if not diff['stats_enabled']:
        content='<div class="notice">未填写 contrasts.tsv；本次仅生成组别/QC描述，不进行差异表达或 DTU 推断。</div>'
    else:
        content+='<p>每个比较按 contrasts.tsv 执行；基因和转录本差异表达使用 edgeR robust quasi-likelihood GLM，DTU 使用 DRIMSeq 并以 stageR 做分阶段校正。显著性阈值为 FDR 0.05。DTU 的 logit_change 表示转录本使用比例的模型系数，正负号表示相对参考组的变化方向。</p>'
        content+=table(comparisons,[('contrast_id','比较 ID'),('test_group','测试组'),('reference_group','参考组'),('samples','纳入样本')],'contrast-summary',True)
        low_power=[c['contrast_id'] for c in diff['contrasts'] if c.get('low_power_warning')]
        if low_power:content+='<div class="notice">以下比较至少一组只有 2 个重复，统计功效可能有限：'+esc(', '.join(low_power))+'。</div>'
        if de_gene:
            top=sorted(de_gene,key=lambda r:pvalue_key(r,'FDR'))[:200]
            content+='<h3>基因差异表达（FDR 最小的前 200 项）</h3>'+table(top,[('contrast_id','比较'),('feature_id','基因 ID'),('logFC','log2 fold change'),('logCPM','logCPM'),('PValue','P 值'),('FDR','BH FDR')],'de-gene',True)
        else:content+='<p>没有通过表达过滤的基因差异检验结果。</p>'
        if de_tx:
            top=sorted(de_tx,key=lambda r:pvalue_key(r,'FDR'))[:200]
            content+='<h3>转录本差异表达（FDR 最小的前 200 项）</h3>'+table(top,[('contrast_id','比较'),('feature_id','转录本 ID'),('logFC','log2 fold change'),('logCPM','logCPM'),('PValue','P 值'),('FDR','BH FDR')],'de-transcript',True)
        else:content+='<p>没有通过表达过滤的转录本差异检验结果。</p>'
        if dtu_gene:
            top=sorted(dtu_gene,key=lambda r:pvalue_key(r,'stageR_gene'))[:200]
            content+='<h3>差异转录本使用基因（stageR 筛查校正值优先）</h3>'+table(top,[('contrast_id','比较'),('gene_id','基因 ID'),('pvalue','原始 P 值'),('stageR_gene','stageR 基因筛查校正值'),('stageR_significant','stageR 基因筛查显著')],'dtu-gene',True)
        else:content+='<p>没有通过表达过滤、可进行 DTU 检验的基因。</p>'
        if dtu_tx:
            top=sorted(dtu_tx,key=lambda r:pvalue_key(r,'stageR_transcript'))[:200]
            columns=[('contrast_id','比较'),('gene_id','基因 ID'),('feature_id','转录本 ID'),('logit_change','DTU logit change'),('pvalue','原始 P 值'),('stageR_gene','stageR 基因筛查校正值'),('stageR_transcript','stageR 转录本分阶段校正值'),('stageR_significant','stageR 显著')]
            content+='<h3>转录本 DTU 结果（stageR 分阶段校正值优先）</h3>'+table(top,columns,'dtu-transcript',True)
        else:content+='<p>没有可报告的转录本 DTU 结果。</p>'
    sections.append(('differential','5 差异表达与 DTU',content))
    evaluable=[r for r in psi if r['psi']!='NA'];content=f'<p>共输出 {ev["psi_rows"]:,} 条样本–事件记录，其中 {ev["evaluable_psi_rows"]:,} 条达到支持阈值。PSI 由统一目录及 IsoQuant TPM 计算。</p>'
    if evaluable:
        fig,ax=plt.subplots(figsize=(7.4,4.4))
        for i,s in enumerate(samples):
            vals=[float(r['psi']) for r in evaluable if r['sample_id']==s]
            if vals:ax.hist(vals,bins=np.linspace(0,1,21),histtype='step',linewidth=1.8,color=sample_color(s),label=s)
        ax.set(xlabel='PSI',ylabel='Evaluable events');ax.legend(frameon=False)
        content+=figure(fd,'psi_distribution',fig,'可评估剪接事件 PSI 分布','事件相关转录本合计至少 10 reads 才展示 PSI；PSI 属于表达推算结果，事件类型计数可能互相重叠。')
        content+='<h3>支持最高的前 200 条可评估样本–事件记录</h3>'+table(sorted(evaluable,key=lambda r:float(r['support_reads']),reverse=True)[:200],[('sample_id','样本'),('gene_id','基因'),('event_type','事件类型'),('event_id','事件'),('support_reads','支持 reads'),('psi','PSI')],'psi-top',True)
    else:content+='<div class="notice">本次没有达到支持要求的可评估剪接事件；空结果已保留，不解释为无可变剪接。</div>'
    sections.append(('splicing','6 可变剪接',content))
    content=f'<p>本次获得 {ev["terminal_candidates"]:,} 条样本级末端候选记录。按基因、链方向与末端类型聚合，默认窗口 25 bp，至少 5 条不同 read ID 支持。</p>'
    if ends:
        fig,ax=plt.subplots(figsize=(7.4,4.4));pos=np.arange(len(samples));c=Counter((r['sample_id'],r['end_type']) for r in ends)
        for i,kind in enumerate(['5prime','3prime']):ax.bar(pos+(i-.5)*.35,[c[(s,kind)] for s in samples],width=.35,color=COLORS[i+2],label=kind)
        ax.set_xticks(pos,samples);ax.set(ylabel='Candidate terminal clusters');ax.legend(frameon=False)
        content+=figure(fd,'terminal_counts',fig,'末端候选数量','5′为读段比对起点候选，3′另要求 poly(A/T) 支持；同一基因内共享末端的异构体不会被重复当成不同位点。')
        content+='<h3>支持最高的前 200 条末端候选</h3>'+table(sorted(ends,key=lambda r:int(r['support_reads']),reverse=True)[:200],[('sample_id','样本'),('gene_id','基因'),('chrom','染色体'),('strand','链'),('end_type','末端'),('position_1based','位置 1-based'),('support_reads','支持 reads')],'ends-top',True)
    else:content+='<div class="notice">本次无达到默认支持要求的末端候选。</div>'
    content+='<div class="notice">5′比对起点不等于已验证的转录起始位点。3′末端候选不等于已验证的 APA；缺少 CAGE、独立末端测序或实验验证。样本间候选数量受深度和读段完整性影响。</div>'
    models=model_info(q/'filtered.gtf');gene_support=Counter()
    for r in passing:gene_support[r['gene_id']]+=sum(float(r[s]) for s in samples)
    for gi,(gene,_) in enumerate(gene_support.most_common(2)):
        candidates=sorted([r for r in passing if r['gene_id']==gene],key=lambda r:sum(float(r[s]) for s in samples),reverse=True)[:6]
        fig,ax=plt.subplots(figsize=(10,3.5))
        for y,r in enumerate(candidates):
            m=models[r['transcript_id']];ex=m['exons']
            if not ex:continue
            ax.plot([ex[0][0],ex[-1][1]],[y,y],color='#89939e',lw=1)
            for start,end in ex:ax.add_patch(Rectangle((start,y-.18),end-start,.36,color=COLORS[3]))
        ax.set_yticks(range(len(candidates)),[r['transcript_id'] for r in candidates]);ax.set(ylim=(-.6,len(candidates)-.4),xlabel='Genomic coordinate (0-based; exon end exclusive)');ax.ticklabel_format(axis='x',style='plain',useOffset=False)
        content+=figure(fd,'structure_track_'+str(gi+1),fig,'高支持基因结构示例：'+gene,'最多展示该基因支持最高的 6 个过滤后模型；坐标为 GRCm39 基因组坐标，按基因组正向排列，链方向见目录。用于查看结构，不构成功能验证。')
    sections.append(('ends','7 末端候选与结构示例',content))
    versions=[]
    for package in ['snakemake','cutadapt','numpy','pandas','scipy','pysam','matplotlib']:
        try:value=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:value='not available'
        versions.append(dict(software=package,version=value))
    for name,cmd in [('IsoQuant',['isoquant','--version']),('SQANTI3',['sqanti3_qc.py','--version']),('SUPPA2',['suppa.py','--version']),('minimap2',['minimap2','--version']),('samtools',['samtools','--version'])]:
        try:r=subprocess.run(cmd,capture_output=True,text=True,timeout=30);value=(r.stdout or r.stderr).strip().split('\n')[0]
        except (OSError,subprocess.TimeoutExpired):value='not available'
        if name=='SUPPA2':value='2.4 (pinned environment)'
        versions.append(dict(software=name,version=value))
    versions.extend(dict(software=name,version=value) for name,value in diff.get('package_versions',{}).items())
    tsv(d/'software_versions.tsv',['software','version'],versions)
    git={}
    for name,cmd in [('commit',['git','rev-parse','HEAD']),('worktree',['git','status','--short'])]:
        try:
            r=subprocess.run(cmd,capture_output=True,text=True);git[name]=r.stdout.strip() if r.returncode==0 else 'uncommitted initial repository'
        except FileNotFoundError:git[name]=os.environ.get('ISOSEQ_GIT_'+('COMMIT' if name=='commit' else 'WORKTREE'),'unavailable')
    relationship='independent biological replicates' if diff['stats_enabled'] else 'groups/contrasts not supplied; descriptive only'
    js(d/'run_metadata.json',dict(config=x.c,git=git,generated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),sample_relationship=relationship,image=os.environ.get('ISOSEQ_IMAGE','iso-seq-bulk:1.1.0'),CPU_limit=int(os.environ.get('ISOSEQ_CORES','32')),memory_limit_GiB=int(os.environ.get('ISOSEQ_MEMORY_MB','131072'))/1024))
    content='<p>流程：原始质控 → 双端接头与方向处理 → minimap2 剪接比对 → IsoQuant 共享建模 → SQANTI3 结构分类 → 表达定量 → 差异表达/DTU → SUPPA2 剪接与末端候选 → 报告。</p><ul><li>保留的 RNA 方向由建库接头及 poly(A/T) 结构支持，不将文件名中的 ISO-seq 自动认定为 PacBio HiFi。</li><li>未确认所有 reads 的 5′完整性，因此未强制 IsoQuant full-length 假设。</li><li>基因与转录本差异表达按 edgeR quasi-likelihood 负二项模型计算；DTU 按 DRIMSeq Dirichlet-multinomial 模型计算，并以 stageR 进行分阶段 FDR 控制。batch 至少有两个水平且在比较内可估计时作为固定效应。</li><li>统计设计假定样本为独立生物学重复，不支持配对或重复测量设计。每组至少 2 个重复；仅 2 个重复时结果需谨慎解释。</li><li>SUPPA2 PSI 和末端候选仍为描述性结果，不进行差异剪接事件检验；短链、内引发、比对、PCR 扩增及参考过滤都可能影响异构体与末端候选。</li><li>本实验为常规 bulk 长读长分析，不估计新生 RNA、合成/降解速率或半衰期。</li></ul>'
    content+='<details><summary>软件版本与运行信息</summary>'+table(versions,[('software','软件'),('version','版本')],'versions')+'<p>输入与参考只读挂载；Snakemake 管理阶段依赖和断点恢复。日志在 log/，中间文件在 tmp/，具体参数及运行版本保留在交付目录。</p></details>'
    sections.append(('methods','8 方法与解读范围',content))
    status=dict(status='analysis_complete',samples=samples,required_modules=['raw_qc','trim','alignment','read_flow_sankey','isoquant','sqanti3','quantification','differential','events'],statistics=diff['status'],optional_absence=dict(differential='no contrasts requested' if not diff['stats_enabled'] else None,splicing='no evaluable events' if not evaluable else None,terminal='no supported candidates' if not ends else None))
    js(d/'workflow_status.json',status)
    files=[]
    for sample in samples:
        for stage in ['01_raw_qc','02_trim','03_alignment','06_quantification']:
            files.extend(sorted((x.o/sample/stage).glob('*')))
    for stage in ['00_preflight','04_catalog','05_structure_qc','06_quantification','07_events','08_differential']:
        files.extend(p for p in sorted((x.p/stage).rglob('*')) if p.is_file() and p.suffix.lower() not in ['.db','.pickle','.pkl','.html'] and '/aux/' not in str(p))
    files.extend(flow_targets)
    files.extend(sorted(fd.glob('*')));files.extend([d/'software_versions.tsv',d/'run_metadata.json',d/'workflow_status.json'])
    manifest=[dict(category='analysis',source=str(p.resolve()),target=str(p.resolve()),status='up_to_date',description=p.name) for p in files if p.is_file() and p.stat().st_size]
    manifest.append(dict(category='final_report',source=str((d/'report.html').resolve()),target=str((d/'report.html').resolve()),status='copied',description='最终中文离线报告'))
    tsv(d/'delivery_manifest.tsv',['category','source','target','status','description'],manifest)
    nav=''.join(f'<a href="#{id}">{esc(title)}</a>' for id,title,c in sections)
    body=''.join(f'<section class="report-section" id="{id}"><h2>{esc(title)}</h2>{c}</section>' for id,title,c in sections)
    plotly=(Path(__file__).resolve().parents[1]/'assets/plotly-3.5.1.min.js').read_text(encoding='utf-8')
    output=f'<!DOCTYPE html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{esc(x.c["report_title"])}</title><style>{CSS}</style><script>{plotly}</script></head><body><nav class="sidebar"><strong>bulk ISO-seq</strong>{nav}</nav><main>{body}</main><footer class="footer">{esc(x.c.get("report_footer","ISO-seq · 中文离线分析交付"))}</footer><div id="modal" hidden role="dialog" aria-label="放大图形"><button id="close-modal">关闭 / Esc</button><img alt="放大图形"></div><script>{JS}</script></body></html>'
    (d/'report.html').write_text(output,encoding='utf-8')
    print('REPORT',d/'report.html',flush=True)

if __name__=='__main__':main()
