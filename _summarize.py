
from kvarq.genes import load_testsuite
from kvarq.analyse import Analyser

# rows : genomes
# columns : SNPs
#
# table1 : mean coverage
# table2 : sd/mean
# table3 : [AX][M]
#   - D if derived has majority
#   - X if non-ancestral-non-derived has majority
#   - M if most prevalent < 90%

import json, csv, sys, os.path

analyser = Analyser()

coll_path = os.path.join(os.path.dirname(__file__), 'coll14.py')
name, testsuite = load_testsuite(coll_path)
tests = testsuite.tests
testsuites = {name: testsuite}

columns = ['filename']
columns += [str(test) for test in tests]

means = {}
sds = {}
types = {}

for fname in sys.argv[1:]:

    sys.stderr.write('reading %s...' % fname)
    d = json.load(file(fname))
    analyser.decode(testsuites, d)

    means[fname] = {'filename':fname}
    sds[fname] = {'filename':fname}
    types[fname] = {'filename':fname}
    #data[fname]['filesize'] = sum(d['info']['size'])
    #data[fname]['scantime'] = int(d['info']['scantime'])

    total = 0
    for test in tests:

        coverage = analyser[test]

        mean = int(coverage.mean())
        total += mean
        means[fname][str(test)] = mean
        sds[fname][str(test)] = '%.2f' % (coverage.std() / max(1, coverage.mean()))

        most = 0
        base = None
        depth = 0
        for b, n in coverage.bases_at(coverage.start).items():
            if n > most:
                most = n
                base = b
            depth += n

        t = ''
        if base == test.template.base:
            t += 'D'
        elif base != test.template.orig:
            t += 'X'
        if most < 0.9 * depth:
            t += 'M'
        types[fname][str(test)] = t

    average = float(total) / len(tests)
    for test in tests:
        means[fname][str(test)] = '%.2f' % (
                means[fname][str(test)] / average
            )

    sys.stderr.write('done\n')

out = csv.writer(sys.stdout)
out.writerow(columns)

out.writerow(['MEAN COVERAGE (fraction of average)'])
for fname in sys.argv[1:]:
    out.writerow([means[fname].get(k) for k in columns])
out.writerow(['SD/MEAN'])
for fname in sys.argv[1:]:
    out.writerow([sds[fname].get(k) for k in columns])
out.writerow(['Derived X-not-derived-not-ancestral Mixed'])
for fname in sys.argv[1:]:
    out.writerow([types[fname].get(k) for k in columns])


