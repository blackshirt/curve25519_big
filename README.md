curve25519
__________

what's this repository for
**************************

This repository contains two files implements `x25519` function
for implementing curve25519 from RFC 7748.
1. in the form 'x25519.v' file, it's backed by standard 'math.big' library of standard vlib.
but, the test included show it's buggy and run terribely slow, hogs your ram and cpu.
2. in the form 'x25519_gmp.v' file and backed by 'v_gmp'. the test run successfully
