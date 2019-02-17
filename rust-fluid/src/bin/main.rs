//use rust_fluid::integration::euler;
//use rust_fluid::integration::bogacki_shampine;
//use rust_fluid::integration::runge_kutta_4;
use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::advection;
use rust_fluid::interpolation;
use rust_fluid::integration;
use rust_fluid::linear_solvers;


fn main() {
//    let mut euler_y = 1.0;
//    let mut bogas_y = 1.0;
//    let mut runge_y = 1.0;
//
//    let dt = 0.02;
//    let f = |t: f64, y: f64| y;
//
//    for i in 1..=50 {
//        let t = i as f64 * dt;
//
//        euler_y = euler(t, euler_y, &f, dt);
//        bogas_y = bogacki_shampine(t, bogas_y, &f, dt);
//        runge_y = runge_kutta_4(t, runge_y, &f, dt);
//    }
//
//    println!("exact: {}", std::f64::consts::E);
//    println!("euler: {}", euler_y);
//    println!("bogas: {}", bogas_y);
//    println!("runge: {}", runge_y);

    let mut buffer = vec![0u8; 512 * 1024 * 4];

    let mut solver = FluidSolver::new(512, 1024, 0.01, 1.0 / 16.0, 1.0)
        .interpolation(interpolation::bilinear_interpolate)
        .integration(integration::bogacki_shampine)
        .linear_solver(linear_solvers::gauss_siedel)
        .advection(advection::semi_lagrangian);

    for i in 0..32 {
        for j in 0..32 {
            *solver.y_velocity_src.at_mut(220 + i, 160 + j) = 150.0;
            *solver.density_src.at_mut(220 + i, 160 + j) = 2.0;
        }
    }


    for iteration in 0..400 {


        solver.solve();

        for i in 0..32 {
            for j in 0..32 {
                *solver.y_velocity_src.at_mut(220 + i, 160 + j) = 100.0;
                *solver.density_src.at_mut(220 + i, 160 + j) = 2.0;
            }
        }

        solver.to_image(&mut buffer);
        lodepng::encode32_file(format!("output/{}.png", iteration), &buffer, 1024, 512).unwrap();



    }

}
