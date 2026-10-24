#!/usr/bin/env python3
"""End-to-end tests using a small independently written multi-resolution Hi-C.

Only Python's standard library is required. No downloaded biological dataset or
mock solver/reader is used: all production executables read these binary maps.
"""
import math
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import zlib


def pack(fmt, *values):
    return struct.pack('<' + fmt, *values)


def cstring(text):
    return text.encode() + b'\0'


def fixture(path, sparse=False, genome='custom', chromosome_lengths=None, resolutions=None):
    resolutions = resolutions or [500, 1000, 7000, 50000, 100000]
    chroms = chromosome_lengths or [('1', 3_200_123), ('Y', 3_200_000), ('MT', 3_200_123), ('Un_test', 3_200_123)]
    if sparse:
        chroms.append(('Tiny', 250))
    data = bytearray(b'HIC\0' + pack('iQ', 9, 0) + cstring(genome))
    data += pack('QQii', 0, 0, 0, len(chroms))
    for name, length in chroms:
        data += cstring(name) + pack('Q', length)
    data += pack('i', len(resolutions)) + pack(f'{len(resolutions)}i', *resolutions) + pack('i', 0)
    matrices = []
    for first in range(len(chroms)):
        for second in range(first, len(chroms)):
            blocks = []
            for resolution in resolutions:
                # Spread a 32x32 contact matrix across the entire chromosome.
                # Distinct nonzero bins occur in both compartments at every res.
                bins = [(i * 3_200_000 // 32) // resolution for i in range(32)]
                records = []
                for y in range(32):
                    for x in range(32 if first != second else y + 1):
                        value = 2 + ((x*13 + y*7 + x*y*3) % 19)/10
                        if (x < 16) == (y < 16):
                            value += 8
                        records.append((bins[x], bins[y], value))
                if chroms[first][0] == 'Tiny' or chroms[second][0] == 'Tiny':
                    records = []
                raw = bytearray(pack('iii4Bh', len(records), 0, 0, 1, 0, 0, 1, 32 if records else 0))
                for y in range(32 if records else 0):
                    row = [(x, value) for x, yy, value in records if yy == bins[y]]
                    raw += pack('hh', bins[y], len(row))
                    for x, value in row:
                        raw += pack('hf', x, value)
                compressed = zlib.compress(raw)
                position = len(data)
                data += compressed
                blocks.append((resolution, position, len(compressed), sum(v for _, _, v in records)))
            position = len(data)
            data += pack('iii', first, second, len(resolutions))
            for resolution, block_at, size, total in blocks:
                data += cstring('BP') + pack('i4f4i', 0, total, 0, 0, 0, resolution, 10000, 1, 1)
                data += pack('iqi', 0, block_at, size)
            matrices.append((f'{first}_{second}', position, len(data)-position))
    master = len(data)
    struct.pack_into('<Q', data, 8, master)
    footer = bytearray(pack('i', len(matrices)))
    for key, position, size in matrices:
        footer += cstring(key) + pack('qi', position, size)
    # Observed/NONE queries do not need expectation or normalization vectors.
    data += pack('Q', len(footer)) + footer + pack('iii', 0, 0, 0)
    path.write_bytes(data)
    return chroms


def run(command, ok=True):
    result = subprocess.run(list(map(str, command)), text=True, capture_output=True, timeout=40)
    if (result.returncode == 0) != ok:
        raise AssertionError(f'{command}\nstatus={result.returncode}\n{result.stdout}\n{result.stderr}')
    return result


def wig(path):
    sections = {}
    current = None
    for line in path.read_text().splitlines():
        if line.startswith('track'):
            continue
        if line.startswith('fixedStep'):
            fields = dict(field.split('=', 1) for field in line.split()[1:])
            assert int(fields['step']) > 0 and int(fields['span']) > 0, line
            current = fields['chrom']
            sections.setdefault(current, [])
        else:
            value = float(line)
            assert math.isfinite(value)
            sections[current].append(value)
    return sections


def oriented(values, resolution):
    first = [x for i, x in enumerate(values) if i*resolution < 1_600_000 and x != 0]
    second = [x for i, x in enumerate(values) if i*resolution >= 1_600_000 and x != 0]
    assert first and second, 'nonzero calls required in both compartments'
    assert sum(first)/len(first) > sum(second)/len(second), 'negative activity correlation'


def main():
    build = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix='possum-cli-') as directory:
        root = Path(directory)
        hic = root/'fixture.hic'
        chroms = fixture(hic)
        track = root/'activity.bedgraph'
        track.write_text(''.join(f'chr{c} 0 1600000 1\nchr{c} 1600000 {length} 0\n'
                                 for c, length in chroms))
        for resolution in (500, 1000, 7000, 100000):
            for chromosome in ('chr1', 'chrY', 'chrM', 'chrUn_test'):
                output = root/'power.txt'
                run([build/'ev.exe', '-v', 0, '-o', '-p', track, hic, chromosome, output, resolution])
                values = [float(x) for x in output.read_text().splitlines()]
                oriented([x if math.isfinite(x) else 0 for x in values], resolution)
                base = root/'lanczos'
                run([build/'LanChr.exe', '-v', 0, '-o', '-n', 'NONE', '-p', track,
                     hic, chromosome, base, resolution])
                sections = wig(base.with_suffix('.Ev1.wig'))
                oriented(next(iter(sections.values())), resolution)
            power = root/'genome_power.wig'
            run([build/'GWevIntra.exe', '-n', 'NONE', '-v', 0, '-o', 'observed', '-p', track, hic, power, resolution])
            sections = wig(power)
            assert set(sections) == {'chr1', 'chrY', 'chrM', 'chrUn_test'}
            for c, length in chroms:
                values = sections['chrM' if c == 'MT' else 'chr'+c]
                assert len(values) == math.ceil(length/resolution)
                oriented(values, resolution)
            base = root/'genome_lanczos'
            run([build/'Lan.exe', '-v', 0, '-o', '-n', 'NONE', '-p', track, '-j', 2, hic, base, resolution])
            sections = wig(base.with_suffix('.wig'))
            assert set(sections) == {'chr1', 'chrY', 'chrM', 'chrUn_test'}
            for values in sections.values():
                oriented(values, resolution)
        # Driver alias selection, 50kb reference orientation, NV=10 output,
        # eigenvalue diagnostics and rescue sign alignment all run in production.
        base = root/'rescue'
        run([build/'Lan.exe', '--rescue', '--chromosomes', 'chr1,chrY,chrM', '-o', '-n', 'NONE',
             '-p', track, '-v', 0, '-j', 2, hic, base, 500])
        assert set(wig(base.with_suffix('.wig'))) == {'chr1', 'chrY', 'chrM'}
        for values in wig(base.with_suffix('.wig')).values():
            oriented(values, 500)
        assert len(base.with_suffix('.rescue.tsv').read_text().splitlines()) == 4
        assert not base.with_suffix('.Ev11.wig').exists()
        # Trans pipeline must use a single global sign and publish only nv EVs.
        base = root/'trans'
        run([build/'GWev.exe', '-v', 0, '-p', track, hic, base, 1000, 2])
        assert set(wig(root/'trans_Ev1.wig')) == {'chr1', 'chrY', 'chrM', 'chrUn_test'}
        assert not (root/'trans_Ev3.wig').exists()
        sparse_hic = root/'sparse.hic'
        fixture(sparse_hic, sparse=True)
        sparse_base = root/'sparse'
        result = run([build/'Lan.exe', '--rescue', '--chromosomes', 'Tiny', '-o', '-n', 'NONE',
                      '-v', 0, sparse_hic, sparse_base, 500])
        assert wig(sparse_base.with_suffix('.wig')) == {'chrTiny': [0.0]}
        assert 'no_high_resolution_eigenvector' in sparse_base.with_suffix('.rescue.tsv').read_text()
        assert 'zero-filled' in result.stderr
        for executable, args in [('ev.exe', [hic, '1', root/'bad', 0]),
                                  ('LanChr.exe', [hic, '1', root/'bad', 0]),
                                  ('GWevIntra.exe', [hic, root/'bad', 0]),
                                  ('GWev.exe', [hic, root/'bad', 0]),
                                  ('Lan.exe', [hic, root/'bad', 0])]:
            run([build/executable, *args], ok=False)
        run([build/'Lan.exe', hic, root/'bad', '500bp'], ok=False)
        run([build/'Lan.exe', '-p', root/'missing.bedgraph', hic, root/'bad', 500], ok=False)
        dump = root/'dump.txt'
        dump.write_text('0 0 3\n0 500 3\n0 1000 1\n0 1500 1\n500 500 3\n500 1000 1\n'
                        '500 1500 1\n1000 1000 3\n1000 1500 3\n1500 1500 3\n')
        run([build/'mainEigen.exe', '-v', 0, '-m', 2, dump, root/'ascii.txt', 500])
        assert len((root/'ascii.txt').read_text().splitlines()) == 4
        dump.write_text('1000 0 1\n0 0 2\n1000 1000 3\n')
        run([build/'mainEigen.exe', '-v', 0, dump, root/'ascii.txt', 500])
        assert len((root/'ascii.txt').read_text().splitlines()) == 3  # max of both coordinates
        dump.write_text('0 0 1\n0 500')
        run([build/'mainEigen.exe', '-v', 0, dump, root/'ascii.txt', 500], ok=False)
        dump.write_text('0 0 1\ninvalid\n')
        run([build/'mainEigen.exe', '-v', 0, dump, root/'ascii.txt', 500], ok=False)
        print('Production CLI tests passed: 500bp/1kb/7kb/100kb, all chromosome aliases, both methods, rescue, trans, ASCII and invalid inputs')


if __name__ == '__main__':
    main()
