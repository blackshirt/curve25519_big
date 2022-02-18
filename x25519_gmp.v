module x25519

import gmp

const (
	g_zero            = gmp.from_u64(0)
	g_one             = gmp.from_u64(1)
	g_two             = gmp.from_u64(2)

	// curve modulo prime for x25519, 2**255 - 19
	gmp_modprime_25519 = gmp.from_str('57896044618658097711785492504343953926634992332820282019728792003956564819949') or {
		panic(err)
	}
	// The constant a24 for x25519
	gmp_a24const_25519 = gmp.from_u64(121665)
		// curve modulo prime for x448 , 2^448 - 2^224 - 1
		/*
		cv_modprime_448   = gmp.from_str('726838724295606890549323807888004534353641360687318060281490199180612328166730772686396383698676545930088884461843637361053498018365439') or {
		panic(err)
	}*/ // 2^448 - 2^224 - 1
		// The constant a24 for x448
		// cv_a24const_448   = gmp.from_u64(39081)
)

pub fn gmp_x25519(mut k []byte, mut u []byte) ?[]byte {
	// first decode k and u and then perform scalar multiply
	// based on formulas from [montgomery].  All calculations are performed modulo p,
	// and then encode to array of bytes
	if k.len != 32 {
		return error('k should 32 for x25519_gmp, but curve has $k.len')
	}
	//_ = k[31]
	mut key := gmp_decode_scalar(mut k) ?
	mut ucord := gmp_decode_x_coordinate(mut u) ?

	mut res := gmp_scalar_multiply(key, ucord) ?

	return gmp_encode_x_coordinate(mut res)
}

fn gmp_decode_scalar(mut k []byte) ?gmp.Bigint {
	if k.len != x25519.scalar_size {
		return error('Not allowed of k $k.len length')
	}

	// For X25519, in order to decode 32 random bytes as an integer scalar,
	// set the three least significant bits of the first byte and the most significant bit
	// of the last to g_zero, set the second most significant bit of the last
	// byte to 1 and, finally, decode as little-endian
	_ = k[31]
	k[0] &= 248
	k[31] &= 127
	k[31] |= 64

	return gmp_decode_little_endian(k)
}

fn gmp_decode_little_endian(b []byte) ?gmp.Bigint {
	// check if b.len matching with bits
	if b.len != x25519.scalar_size {
		return error('Not allowed b $b.len length')
	}

	//_ = b[55]
	mut sum := x25519.g_zero
	for i in 0 .. (x25519.bits_size + 7) / 8 {
		// mut val := gmp.new()
		// left shift
		// reflect changes of v_gmp lib
		val := gmp.mul_2exp(gmp.from_u64(b[i]), u64(8 * i))
		//val := gmp.from_u64(b[i]).lshift(u32(8 * i))
		sum = sum + val
	}
	return sum
}

fn gmp_encode_x_coordinate(mut u gmp.Bigint) []byte {
	mut arr := []byte{len: (x25519.bits_size + 7) / 8}
	// u = u % gmp_modprime_25519
	u = gmp.mod(u, gmp_modprime_25519)
	//u = gmp.mod(u, x25519.gmp_modprime_25519)
	//_, u = u.div_mod(x25519.gmp_modprime_25519)
	// placed this allocation out of for loop and call gmp.clear() to release it
	mut val := gmp.new()
	for i in 0 .. arr.len {
		// this do right shifting
		val = gmp.tdiv_q_2exp(u, u64(8 * i))
		//val = u.rshift(u32(8 * i))
		// mut d := gmp.new()
		val = gmp.and(val, gmp.from_u64(0xff))
		// arr[i] = byte(d.u64())

		//val = gmp.and(val, gmp.from_u64(0xff))
		arr[i] = byte(val.u64())
	}
	// gmp.clear(mut val)

	return arr
}

fn gmp_decode_x_coordinate(mut u []byte) ?gmp.Bigint {
	// When receiving such an array, implementations of X25519
	//(but not X448) MUST mask the most significant bit in the final byte
	if x25519.bits_size % 8 != 0 {
		u[u.len - 1] &= (1 << (x25519.bits_size % 8)) - 1
	}
	return gmp_decode_little_endian(u)
}

fn gmp_scalar_multiply(k gmp.Bigint, u gmp.Bigint) ?gmp.Bigint {
	mut x_1 := u
	mut x_2 := x25519.g_one
	mut z_2 := x25519.g_zero
	mut x_3 := u
	mut z_3 := x25519.g_one
	mut swap := gmp.new()

	// allocate bigint vars needed for operation outside the loop
	
	mut k_t := gmp.new()
	mut val := gmp.new()
	mut aa := gmp.new()
	mut bb := gmp.new()
	mut xx := gmp.new()
	for t := x25519.bits_size - 1; t >= 0; t-- {
		val = gmp.tdiv_q_2exp(k, u64(t)) // right shift

		k_t = gmp.and(val, x25519.g_one)
		val = gmp.xor(swap, k_t)
		
		swap = val

		x_2, x_3 = gmp_cswap(swap, x_2, x_3)
		z_2, z_3 = gmp_cswap(swap, z_2, z_3)
		swap = k_t

		mut a := (x_2 + z_2) % x25519.gmp_modprime_25519

		aa = gmp.powm_sec(a, g_two, gmp_modprime_25519)
		

		mut b := (x_2 - z_2) % x25519.gmp_modprime_25519

		bb = gmp.powm_sec(b, g_two, gmp_modprime_25519)
		
		mut e := (aa - bb) % x25519.gmp_modprime_25519
		mut c := (x_3 + z_3) % x25519.gmp_modprime_25519
		mut d := (x_3 - z_3) % x25519.gmp_modprime_25519

		mut da := (d * a) % x25519.gmp_modprime_25519
		mut cb := (c * b) % x25519.gmp_modprime_25519

		x_3 = gmp.powm_sec((da + cb) % gmp_modprime_25519, g_two, gmp_modprime_25519)
		

		xx = gmp.powm_sec((da - cb) % gmp_modprime_25519, g_two, gmp_modprime_25519)
		
		z_3 = (x_1 * xx) % x25519.gmp_modprime_25519
		x_2 = (aa * bb) % x25519.gmp_modprime_25519
		z_2 = e * ((aa + (x25519.gmp_a24const_25519 * e) % x25519.gmp_modprime_25519) % x25519.gmp_modprime_25519)
	}
	x_2, x_3 = gmp_cswap(swap, x_2, x_3)
	z_2, z_3 = gmp_cswap(swap, z_2, z_3)

	val = gmp.powm_sec(z_2, gmp_modprime_25519 - g_two, gmp_modprime_25519)

	res := (x_2 * val) % x25519.gmp_modprime_25519

	return res
}

fn gmp_cswap(swap gmp.Bigint, x_2 gmp.Bigint, x_3 gmp.Bigint) (gmp.Bigint, gmp.Bigint) {
	mask := g_zero - swap
	
	vv := gmp.xor(x_2, x_3)
	dummy := gmp.and(mask, vv)
	
	r1 := gmp.xor(x_2, dummy)
	r2 := gmp.xor(x_3, dummy)

	return r1, r2
}
