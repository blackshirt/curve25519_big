module x25519

import encoding.hex
import math.big

struct Xdata {
	k   string
	u   string
	out string
}

fn test_curve_x25519_rfc_vector() ? {
	data := [
		Xdata{
			k: 'a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4'
			u: 'e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c'
			out: 'c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552'
		},
		Xdata{
			k: '4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d'
			u: 'e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493'
			out: '95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957'
		},
	]

	for t in data {
		mut scalar := hex.decode(t.k) ?
		mut cord := hex.decode(t.u) ?

		expected := hex.decode(t.out) ?

		result := x25519(mut scalar, mut cord) ?
		assert result == expected
	}
}

struct ScalMul {
	num1   string
	num2   string
	result string
}

/*
fn test_scalar_multiply() ? {
	data := [
		ScalMul{
			num1: '34426434033919594451155107781188821651316167215306631574996226621102155684838'
			num2: '8883857351183929894090759386610649319417338800022198945255395922347792736741'
			result: '30288769775299154322453471711435182976155832952747871841100586939672915661451'
		},
		// from working python implementation
		ScalMul{
			num1: '31029842492115040904895560451863089656472772604678260265531221036453811406496'
			num2: '35156891815674817266734212754503633747128614016119564763269015315466259359304'
			result: '9687847192437982599224399286282328413310936304537917893524998943417340894289'
		},
	]
	for t in data {
		n1 := big.integer_from_string(t.num1) or { panic(err) }
		n2 := big.integer_from_string(t.num2) or { panic(err) }
		exp := big.integer_from_string(t.result) or { panic(err) }
		// cv := new_curve(255) ?
		res := scalar_multiply(n1, n2) ?
		// assert gmp.cmp(res, exp) == 0
		assert res == exp
	}
}
*/
