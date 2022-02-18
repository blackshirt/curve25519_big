module x25519

import math.big

const (
	cv_zero           = big.integer_from_u64(0)
	cv_one            = big.integer_from_u64(1)
	cv_two            = big.integer_from_u64(2)

	// curve modulo prime for x25519, 2**255 - 19
	cv_modprime_25519 = big.integer_from_string('57896044618658097711785492504343953926634992332820282019728792003956564819949') or {
		panic(err)
	}
	// The constant a24 for x25519
	cv_a24const_25519 = big.integer_from_u64(121665)
		// curve modulo prime for x448 , 2^448 - 2^224 - 1
		/*
		cv_modprime_448   = big.integer_from_string('726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018365439') or {
		panic(err)
	}*/ // 2^448 - 2^224 - 1
		// The constant a24 for x448
		// cv_a24const_448   = big.integer_from_u64(39081)
)

const (
	scalar_size = 32
	bits_size   = 255
)

pub fn x25519(mut k []byte, mut u []byte) ?[]byte {
	// first decode k and u and then perform scalar multiply
	// based on formulas from [montgomery].  All calculations are performed modulo p,
	// and then encode to array of bytes
	if k.len != 32 {
		return error('k should 32 for x25519, but curve has $k.len')
	}
	//_ = k[31]
	mut key := decode_scalar(mut k) ?
	mut ucord := decode_x_coordinate(mut u) ?

	mut res := scalar_multiply(key, ucord) ?

	return encode_x_coordinate(mut res)
}

fn decode_scalar(mut k []byte) ?big.Integer {
	if k.len != x25519.scalar_size {
		return error('Not allowed of k $k.len length')
	}

	// For X25519, in order to decode 32 random bytes as an integer scalar,
	// set the three least significant bits of the first byte and the most significant bit
	// of the last to cv_zero, set the second most significant bit of the last
	// byte to 1 and, finally, decode as little-endian
	_ = k[31]
	k[0] &= 248
	k[31] &= 127
	k[31] |= 64

	return decode_little_endian(k)
}

fn decode_little_endian(b []byte) ?big.Integer {
	// check if b.len matching with bits
	if b.len != x25519.scalar_size {
		return error('Not allowed b $b.len length')
	}

	//_ = b[55]
	mut sum := x25519.cv_zero
	for i in 0 .. (x25519.bits_size + 7) / 8 {
		// mut val := gmp.new()
		// left shift
		// reflect changes of v_gmp lib
		// val := big.mul_2exp(big.integer_from_u64(b[i]), u64(8 * i))
		val := big.integer_from_u64(b[i]).lshift(u32(8 * i))
		sum = sum + val
	}
	return sum
}

fn encode_x_coordinate(mut u big.Integer) []byte {
	mut arr := []byte{len: (x25519.bits_size + 7) / 8}
	// u = u % cv_modprime_25519
	// gmp.mod(mut u, u, cv_modprime_25519)
	// u = big.mod(u, cv_modprime_25519)
	_, u = u.div_mod(x25519.cv_modprime_25519)
	// placed this allocation out of for loop and call gmp.clear() to release it
	mut val := big.Integer{}
	for i in 0 .. arr.len {
		// this do right shifting
		// val = big.tdiv_q_2exp(u, u64(8 * i))
		val = u.rshift(u32(8 * i))
		// mut d := gmp.new()
		// gmp.and(mut d, val, gmp.integer_from_u64(0xff))
		// arr[i] = byte(d.u64())

		val = val.bitwise_and(big.integer_from_u64(0xff))
		arr[i] = byte(val.int())
	}
	// gmp.clear(mut val)

	return arr
}

fn decode_x_coordinate(mut u []byte) ?big.Integer {
	// When receiving such an array, implementations of X25519
	//(but not X448) MUST mask the most significant bit in the final byte
	if x25519.bits_size % 8 != 0 {
		u[u.len - 1] &= (1 << (x25519.bits_size % 8)) - 1
	}
	return decode_little_endian(u)
}

fn scalar_multiply(k big.Integer, u big.Integer) ?big.Integer {
	if k.len != 32 {
		return error("curve size doesn't allowed, get $k.len")
	}
	mut x_1 := u
	mut x_2 := x25519.cv_one
	mut z_2 := x25519.cv_zero
	mut x_3 := u
	mut z_3 := x25519.cv_one
	mut swap := u64(0)

	// allocate bigint vars needed for operation outside the loop
	// mut k_t := big.Integer{}
	mut k_t := u64(0)
	mut val := big.Integer{}
	mut aa := big.Integer{}
	mut bb := big.Integer{}
	mut xx := big.Integer{}
	for t := x25519.bits_size - 1; t >= 0; t-- {
		// val = big.tdiv_q_2exp(k, u64(t)) // right shift
		k_t = u64((k.rshift(u32(t)).int()) & 1)
		swap ^= k_t
		// k_t = big.and(val, cv_one)

		// val = big.xor(swap, k_t)

		// swap = val

		x_2, x_3 = cswap(swap, mut x_2, mut x_3)
		z_2, z_3 = cswap(swap, mut z_2, mut z_3)
		swap = k_t

		mut a := (x_2 + z_2) % x25519.cv_modprime_25519

		// this is not available on big.Integer
		// aa = big.powm_sec(a, cv_two, cv_modprime_25519)
		aa = a.mod_pow(u32(2), x25519.cv_modprime_25519)

		mut b := (x_2 - z_2) % x25519.cv_modprime_25519

		// bb = big.powm_sec(b, cv_two, cv_modprime_25519)
		bb = b.mod_pow(u32(2), x25519.cv_modprime_25519)

		mut e := (aa - bb) % x25519.cv_modprime_25519
		mut c := (x_3 + z_3) % x25519.cv_modprime_25519
		mut d := (x_3 - z_3) % x25519.cv_modprime_25519

		mut da := (d * a) % x25519.cv_modprime_25519
		mut cb := (c * b) % x25519.cv_modprime_25519

		// x_3 = big.powm_sec((da + cb) % cv_modprime_25519, cv_two, cv_modprime_25519)
		x_3 = (da + cb) % x25519.cv_modprime_25519
		x_3 = x_3.mod_pow(u32(2), x25519.cv_modprime_25519)

		// xx = big.powm_sec((da - cb) % cv_modprime_25519, cv_two, cv_modprime_25519)
		xx = (da - cb) % x25519.cv_modprime_25519
		xx = xx.mod_pow(u32(2), x25519.cv_modprime_25519)

		// z_3 = (x_1 * xx) % cv_modprime_25519
		z_3 = (x_1 * xx) % x25519.cv_modprime_25519
		x_2 = (aa * bb) % x25519.cv_modprime_25519
		z_2 = e * ((aa + (x25519.cv_a24const_25519 * e) % x25519.cv_modprime_25519) % x25519.cv_modprime_25519)
	}
	x_2, x_3 = cswap(swap, mut x_2, mut x_3)
	z_2, z_3 = cswap(swap, mut z_2, mut z_3)

	// val = big.powm_sec(z_2, cv_modprime_25519 - cv_two, cv_modprime_25519)
	val = z_2.big_mod_pow(x25519.cv_modprime_25519 - x25519.cv_two, x25519.cv_modprime_25519)
	// res := (x_2 * zz) % cv_modprime_25519
	res := (x_2 * val) % x25519.cv_modprime_25519

	return res
}

fn mask_big_bits(cond u64) big.Integer {
	// in go, `^` operates on bit mean NOT, flip bit
	// in v, its a ~    bitwise NOT
	s := ~(u64(cond) - 1)
	return big.integer_from_u64(s)
}

fn cswap(swap u64, mut v big.Integer, mut u big.Integer) (big.Integer, big.Integer) {
	m := mask_big_bits(swap)
	mut t := m.bitwise_and(v.bitwise_xor(u))
	v = v.bitwise_xor(t)
	u = u.bitwise_xor(t)
	return u, v
}
