import sys, tempfile, unittest, json, gzip
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from analyze import Context, orient, fastq, end_clusters, matrix, js
from design import validate_design

class AnalysisTest(unittest.TestCase):
    def design_fixture(self, root, groups, batches=None, contrasts=None):
        config = root / 'config'
        config.mkdir(parents=True, exist_ok=True)
        samples_path = config / 'samples.tsv'
        contrasts_path = config / 'contrasts.tsv'
        with samples_path.open('w', newline='') as handle:
            writer = __import__('csv').writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['sample_id', 'fastq', 'group', 'batch'])
            for i, group in enumerate(groups):
                writer.writerow([f'S{i+1}', f'inputs/S{i+1}.fastq.gz', group, (batches or [''] * len(groups))[i]])
        with contrasts_path.open('w', newline='') as handle:
            writer = __import__('csv').writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['contrast_id', 'test_group', 'reference_group'])
            for row in contrasts or []:
                writer.writerow(row)
        return samples_path, contrasts_path

    def test_multi_group_batch_design_and_explicit_contrasts(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            root = Path(td)
            groups = ['Control'] * 3 + ['DrugA'] * 3 + ['DrugB'] * 3
            batches = ['B1', 'B2', 'B1'] * 3
            samples, contrasts = self.design_fixture(root, groups, batches, [
                ['DrugA_vs_Control', 'DrugA', 'Control'],
                ['DrugB_vs_Control', 'DrugB', 'Control'],
            ])
            result = validate_design(samples, contrasts, root, check_fastq=False)
            self.assertTrue(result['stats_enabled'])
            self.assertTrue(result['batch_enabled'])
            self.assertEqual(len(result['contrasts']), 2)
            self.assertEqual(result['group_counts']['DrugB'], 3)

    def test_two_replicates_warn_but_are_estimable(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            root = Path(td)
            samples, contrasts = self.design_fixture(root, ['A', 'A', 'B', 'B'], contrasts=[['B_vs_A', 'B', 'A']])
            result = validate_design(samples, contrasts, root, check_fastq=False)
            self.assertEqual(result['low_power_groups'], ['A', 'B'])

    def test_batch_confounded_with_group_is_rejected(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            root = Path(td)
            samples, contrasts = self.design_fixture(
                root, ['A', 'A', 'B', 'B'], ['B1', 'B1', 'B2', 'B2'], [['B_vs_A', 'B', 'A']]
            )
            with self.assertRaisesRegex(ValueError, 'confounded'):
                validate_design(samples, contrasts, root, check_fastq=False)

    def test_descriptive_mode_accepts_unassigned_groups(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            root = Path(td)
            samples, contrasts = self.design_fixture(root, ['', ''])
            result = validate_design(samples, contrasts, root, check_fastq=False)
            self.assertFalse(result['stats_enabled'])

    def test_orientation_and_quality(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            root=Path(td);(root/'tmp/S').mkdir(parents=True);(root/'outputs/S/02_trim').mkdir(parents=True)
            core='GCCTAGACGTAC'+'T'*15+'ACGT'*30;quality='!'*12+'I'*(len(core)-12)
            with gzip.open(root/'tmp/S/adapters.fastq.gz','wt') as f:
                f.write('@one\n'+core+'\n+\n'+quality+'\n@two rc\n'+core+'\n+\n'+quality+'\n')
            js(root/'outputs/S/02_trim/cutadapt.json',{'read_counts':{'input':3}})
            x=Context.__new__(Context);x.o=root/'outputs';x.t=root/'tmp';x.c={'poly_t_min':12,'poly_t_window':80,'min_insert_length':50}
            orient(x,'S',1)
            reads=list(fastq(root/'tmp/S/insert.fastq.gz'))
            self.assertEqual(len(reads),2);self.assertEqual(reads[0][0],'S|one');self.assertTrue(reads[0][1].endswith('A'*15));self.assertEqual(reads[0][2],'I'*len(reads[0][1]))
            summary=json.loads((root/'outputs/S/02_trim/trim_summary.json').read_text());self.assertEqual(summary['adapter_unassigned_reads'],1)
    def test_cluster_no_transitive_chaining(self):
        self.assertEqual([len(g) for g in end_clusters([(0,'a'),(20,'b'),(40,'c')],25)],[2,1])
    def test_grouped_matrix_comment_header(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[2]/'tmp') as td:
            p=Path(td)/'counts.tsv';p.write_text('#feature_id\tS1\tS2\nt1\t2\t4\n__no_feature\t8\t6\n')
            self.assertEqual(matrix(p,['S1','S2']),{'t1':{'S1':2.,'S2':4.}})

if __name__=='__main__':unittest.main()
