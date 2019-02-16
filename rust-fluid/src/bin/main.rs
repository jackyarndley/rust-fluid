use rust_fluid::integration::euler;
use rust_fluid::integration::bogacki_shampine;
use rust_fluid::integration::runge_kutta_4;

fn main() {
    let mut euler_y = 1.0;
    let mut bogas_y = 1.0;
    let mut runge_y = 1.0;

    let dt = 0.02;
    let f = |t: f64, y: f64| y;

    for i in 1..=50 {
        let t = i as f64 * dt;

        euler_y = euler(t, euler_y, &f, dt);
        bogas_y = bogacki_shampine(t, bogas_y, &f, dt);
        runge_y = runge_kutta_4(t, runge_y, &f, dt);
    }

    println!("exact: {}", std::f64::consts::E);
    println!("euler: {}", euler_y);
    println!("bogas: {}", bogas_y);
    println!("runge: {}", runge_y);
}
