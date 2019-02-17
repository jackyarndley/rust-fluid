use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::advection;
use rust_fluid::interpolation;
use rust_fluid::integration;
use rust_fluid::linear_solvers;


fn main() {
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
