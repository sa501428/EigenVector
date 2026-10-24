#!/usr/bin/env python3
"""Run production solvers against all four real bundles and independently audit their reports."""
import bisect
import csv
import gzip
import math
from pathlib import Path
import statistics
import sys
import tempfile

from test_cli import fixture, run, wig

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'scripts'))
from phasing_reference import count_chromosomes, read_sizes


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0]*len(values)
    i = 0
    while i < len(order):
        j = i+1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        for k in range(i, j):
            result[order[k]] = (i+j-1)/2
        i = j
    return result


def rho(x, y):
    pairs = [(a,b) for a,b in zip(x,y) if math.isfinite(a) and a != 0 and math.isfinite(b)]
    if len(pairs) < 2:
        return math.nan, len(pairs)
    a,b = map(list, zip(*pairs))
    a,b = ranks(a), ranks(b)
    if max(a) == min(a) or max(b) == min(b):
        return math.nan, len(pairs)
    return statistics.correlation(a,b), len(pairs)


def references(assembly, values, resolution, chromosome='chr1'):
    folder = REPO/'resources'/'compartment_phasing'/assembly
    length = read_sizes(folder/'chrom.sizes')[chromosome]
    data = next(data for chrom,_,data in count_chromosomes(folder/'gc.100bp.counts.bin.gz') if chrom == chromosome)
    spans = []
    with gzip.open(folder/'refseq.coverage.bedGraph.gz', 'rt') as stream:
        for line in stream:
            chrom,start,end,value = line.split()
            if chrom == chromosome and value == '1':
                spans.append((int(start),int(end)))
    ends = [b for a,b in spans]
    gc, genes = [math.nan]*len(values), [math.nan]*len(values)
    for i,x in enumerate(values):
        if not math.isfinite(x) or x == 0:
            continue
        start,end = i*resolution,min(length,(i+1)*resolution)
        if resolution % 100 == 0:
            first,last = start//100,(end+99)//100
            g,v = sum(data[2*first:2*last:2]),sum(data[2*first+1:2*last:2])
            if v and v/(end-start) >= 0.5:
                gc[i] = g/v
        covered = 0
        j = bisect.bisect_right(ends,start)
        while j < len(spans) and spans[j][0] < end:
            covered += min(end,spans[j][1])-max(start,spans[j][0])
            j += 1
        genes[i] = covered/(end-start)
    return gc,genes


def audit(path, values, assembly, resolution):
    with open(path) as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    assert len(rows) == 1, rows
    row = rows[0]
    assert row['assembly'] == assembly and int(row['resolution']) == resolution, row
    sign = int(row['applied_sign']) or 1
    gc,genes = references(assembly,values,resolution)
    for name,reference in [('gc',gc),('refseq',genes)]:
        actual,count = rho(values,reference)
        assert int(row[name+'_bins']) == count, (name,row,count)
        recorded = float(row[name+'_spearman'])
        assert (math.isnan(actual) and math.isnan(recorded)) or abs(actual-recorded*sign) < 1e-8, (name,row,actual)
    strong_gc = int(row['gc_bins']) >= 20 and abs(float(row['gc_spearman'])) >= 0.1
    strong_genes = int(row['refseq_bins']) >= 20 and abs(float(row['refseq_spearman'])) >= 0.1
    if strong_gc and strong_genes and float(row['gc_spearman'])*float(row['refseq_spearman']) < 0:
        assert row['status'] == 'conflicting_references' and row['decision_sign'] == '0' and row['applied_sign'] == '1', row
    elif strong_gc or strong_genes:
        expected = 'bundled_gc' if strong_gc else 'bundled_refseq'
        assert row['source'] == expected and row['status'] == 'oriented', row
        assert float(row['gc_spearman' if strong_gc else 'refseq_spearman'])*sign > 0
    else:
        assert row['decision_sign'] == '0' and row['applied_sign'] == '1', row
    return row


def main():
    build = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix='possum-bundled-') as directory:
        root = Path(directory)
        aliases = {'hg19':'project_hg19_custom','hg38':'project_HG38_no_alt','mm9':'mm9_sample','mm10':'build_mm10_custom'}
        for assembly,label in aliases.items():
            length = read_sizes(REPO/'resources'/'compartment_phasing'/assembly/'chrom.sizes')['chr1']
            hic = root/f'{assembly}.hic'
            fixture(hic,genome=label,chromosome_lengths=[('1',length)],resolutions=[7000])
            output = root/f'{assembly}.power.txt'
            run([build/'ev.exe','-v',0,'-o',hic,'1',output,7000])
            values = [float(line) for line in output.read_text().splitlines()]
            audit(str(output)+'.orientation.tsv',values,assembly,7000)
            base = root/f'{assembly}.lanczos'
            run([build/'LanChr.exe','-v',0,'-o','-n','NONE',hic,'1',base,7000])
            audit(str(base)+'.orientation.tsv',wig(Path(str(base)+'.Ev1.wig'))['chr1'],assembly,7000)
            genome = root/f'{assembly}.genome.wig'
            run([build/'GWevIntra.exe','-n','NONE','-v',0,'-o','observed',hic,genome,7000])
            audit(str(genome)+'.orientation.tsv',wig(genome)['chr1'],assembly,7000)
            base = root/f'{assembly}.driver'
            run([build/'Lan.exe','-v',0,'-o','-n','NONE',hic,base,7000])
            audit(str(base)+'.orientation.tsv',wig(Path(str(base)+'.wig'))['chr1'],assembly,7000)
            fallback_hic = root/f'{assembly}.fallback.hic'
            fixture(fallback_hic,genome='custom_reference',chromosome_lengths=[('1',length)],resolutions=[7000])
            run([build/'ev.exe','-v',0,'-o',fallback_hic,'1',output,7000])
            values = [float(line) for line in output.read_text().splitlines()]
            row = audit(str(output)+'.orientation.tsv',values,assembly,7000)
            assert row['assembly_source'] == 'chromosome_sizes'
            print(f'PASS {assembly}: both methods, chromosome/genome entry points, independent GC/RefSeq Spearman audit',flush=True)
        # A true trans calculation must preserve relative chromosome signs.
        sizes = read_sizes(REPO/'resources'/'compartment_phasing'/'hg38'/'chrom.sizes')
        trans_hic = root/'trans.hic'
        fixture(trans_hic,genome='hg38',chromosome_lengths=[('1',sizes['chr1']),('Y',sizes['chrY'])],resolutions=[7000])
        raw_base,phased_base = root/'trans_raw',root/'trans_phased'
        run([build/'GWev.exe','--phasing-method','none','-v',0,trans_hic,raw_base,7000,1])
        run([build/'GWev.exe','-v',0,trans_hic,phased_base,7000,1])
        raw,phased = wig(Path(str(raw_base)+'_Ev1.wig')),wig(Path(str(phased_base)+'_Ev1.wig'))
        with open(str(phased_base)+'.orientation.tsv') as stream:
            rows = list(csv.DictReader(stream,delimiter='\t'))
        assert len(rows) == 1 and rows[0]['chromosome'] == 'genome'
        row = rows[0]
        sign = int(row['applied_sign'])
        all_values,all_gc,all_genes = [],[],[]
        for chrom in phased:
            assert all(abs(a-b*sign) < 1e-9 for a,b in zip(phased[chrom],raw[chrom]))
            gc,genes = references('hg38',phased[chrom],7000,chrom)
            all_values.extend(phased[chrom]); all_gc.extend(gc); all_genes.extend(genes)
        for name,ref in [('gc',all_gc),('refseq',all_genes)]:
            actual,count = rho(all_values,ref)
            assert count == int(row[name+'_bins'])
            recorded = float(row[name+'_spearman'])
            assert (math.isnan(actual) and math.isnan(recorded)) or abs(actual-recorded*sign) < 1e-8
        # Rescue forwards all controls, retains candidate reports and cleans up.
        rescue_hic = root/'rescue.hic'
        fixture(rescue_hic,genome='hg38',chromosome_lengths=[('1',sizes['chr1'])],resolutions=[7000,50000])
        base = root/'rescued'
        run([build/'Lan.exe','--rescue','--phasing-method','refseq','--phasing-min-correlation','0',
             '--phasing-min-bins','2','-v',0,'-o','-n','NONE',rescue_hic,base,7000])
        with open(str(base)+'.orientation.tsv') as stream:
            row = next(csv.DictReader(stream,delimiter='\t'))
        assert row['assembly'] == 'hg38' and row['gc_bins'] == '0' and row['min_bins'] == '2'
        assert len(Path(str(base)+'.candidates.orientation.tsv').read_text().splitlines()) == 11
        assert len(Path(str(base)+'.lowres_50kb.orientation.tsv').read_text().splitlines()) == 2
        assert not list(root.glob('rescued.lan_tmp.*'))
        print('PASS native trans global sign/correlations and rescue option forwarding/report retention',flush=True)
        length = read_sizes(REPO/'resources'/'compartment_phasing'/'hg38'/'chrom.sizes')['chr1']
        hic = root/'default.hic'
        fixture(hic,genome='',chromosome_lengths=[('1',length)],resolutions=[7777])
        output = root/'default.txt'
        run([build/'ev.exe','-v',0,'-o',hic,'1',output,7777])
        values = [float(line) for line in output.read_text().splitlines()]
        row = audit(str(output)+'.orientation.tsv',values,'hg38',7777)
        assert row['gc_bins'] == '0' and row['gc_exact'] == '0'
        assert row['assembly_source'] == 'chromosome_sizes'
        run([build/'ev.exe','--phasing-method','gc','-o',hic,'1',output,7777],ok=False)
        run([build/'ev.exe','--phasing-method','gc','--phasing-allow-approximate-gc','-o',hic,'1',output,7777])
        with open(str(output)+'.orientation.tsv') as stream:
            row = next(csv.DictReader(stream, delimiter='\t'))
        assert row['gc_exact'] == '0' and int(row['gc_bins']) >= 20
        # Missing, unrecognized, ambiguous and contradictory IDs use the full
        # chromosome metadata when it uniquely identifies a supported build.
        fixture(hic,genome='foreign',chromosome_lengths=[('1',length)],resolutions=[7000])
        run([build/'ev.exe','-o',hic,'1',output,7000])
        values = [float(line) for line in output.read_text().splitlines()]
        row = audit(str(output)+'.orientation.tsv',values,'hg38',7000)
        assert row['assembly_source'] == 'chromosome_sizes'
        run([build/'ev.exe','--assembly','hg38','-o',hic,'1',output,7000])
        values = [float(line) for line in output.read_text().splitlines()]
        audit(str(output)+'.orientation.tsv',values,'hg38',7000)
        for label in ('hg19_wrong_label','hg19_to_hg38'):
            fixture(hic,genome=label,chromosome_lengths=[('1',length)],resolutions=[7000])
            run([build/'ev.exe','-o',hic,'1',output,7000])
            values = [float(line) for line in output.read_text().splitlines()]
            row = audit(str(output)+'.orientation.tsv',values,'hg38',7000)
            assert row['assembly_source'] == 'chromosome_sizes'
        run([build/'ev.exe','--assembly','hg19','-o',hic,'1',output,7000],ok=False)
        hg19_chr2 = read_sizes(REPO/'resources'/'compartment_phasing'/'hg19'/'chrom.sizes')['chr2']
        fixture(hic,genome='foreign',chromosome_lengths=[('1',length),('2',hg19_chr2)],resolutions=[7000])
        run([build/'ev.exe','-o',hic,'1',output,7000])
        assert 'unsupported_assembly' in Path(str(output)+'.orientation.tsv').read_text()
        fixture(hic,genome='hg38',resolutions=[7000])
        run([build/'ev.exe','-o',hic,'1',output,7000],ok=False)
        run([build/'ev.exe','--phasing-method','none','-o',hic,'1',output,7000])
        assert 'disabled' in Path(str(output)+'.orientation.tsv').read_text()
        print('PASS suffixed IDs, all four size fallbacks, whole-map validation, overrides, ambiguity, off-grid fallback/opt-in and disabled mode')


if __name__ == '__main__':
    main()
