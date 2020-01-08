pub fn triangle_occupancy(out1: f64, in1: f64, out2: f64) -> f64 {
    0.5 * in1 * in1 / ((out1 - in1) * (out2 - in1))
}

pub fn trapezoid_occupancy(out1: f64, out2: f64, in1: f64, in2: f64) -> f64 {
    0.5 * (-in1 / (out1 - in1) - in2 / (out2 - in2))
}

pub fn occupancy(d11: f64, d12: f64, d21: f64, d22: f64) -> f64 {
    let mut b: u8 = 0;

    let ds = [d11, d12, d22, d21];

    for i in (0..4).rev() {
        b = (b << 1) | if ds[i] < 0.0 {1} else {0};
    }

    return match b {
        0x0 => 0.0,
        0x1 => triangle_occupancy(d21, d11, d12),
        0x2 => triangle_occupancy(d11, d12, d22),
        0x4 => triangle_occupancy(d12, d22, d21),
        0x8 => triangle_occupancy(d22, d21, d11),
        0xE => 1.0 - triangle_occupancy(-d21, -d11, -d12),
        0xD => 1.0 - triangle_occupancy(-d11, -d12, -d22),
        0xB => 1.0 - triangle_occupancy(-d12, -d22, -d21),
        0x7 => 1.0 - triangle_occupancy(-d22, -d21, -d11),
        0x3 => trapezoid_occupancy(d21, d22, d11, d12),
        0x6 => trapezoid_occupancy(d11, d21, d12, d22),
        0x9 => trapezoid_occupancy(d12, d22, d11, d21),
        0xC => trapezoid_occupancy(d11, d12, d21, d22),
        0x5 => triangle_occupancy(d11, d12, d22) + triangle_occupancy(d22, d21, d11),
        0xA => triangle_occupancy(d21, d11, d12) + triangle_occupancy(d12, d22, d21),
        0xF => 1.0,
        _ => 0.0
    }
}