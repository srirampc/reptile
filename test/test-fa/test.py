import sys

fa_file = sys.argv[1]
q_file = sys.argv[2]
fafp = open(fa_file, 'r')
qfp = open(q_file, 'r')
while True:
    fa_line1 = fafp.readline()
    fa_line2 = fafp.readline()
    q_line1 = qfp.readline()
    q_line2 = qfp.readline()
    if not fa_line2: break  # EOF
    if not q_line2: break  # EOF
    print fa_line1.replace(">", "@"),
    print fa_line2,
    print q_line1.replace(">", "+"),
    print q_line2[1:],

fafp.close()
qfp.close()
