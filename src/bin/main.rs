use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::{integration, linear_solvers, advection};
use rust_fluid::util::sparse::Sparse;
use rust_fluid::solid::SolidBody;

extern crate image;

fn main() {
    let width = 128;
    let height = 128;

    let mut buffer = vec![0u8; width * height * 4];

    let mut bodies = Vec::new();
    bodies.push(SolidBody::new_box(0.5, 0.6, 0.7, 0.1, std::f64::consts::FRAC_PI_4, 0.0, 0.0, 0.0));

    let mut solver = FluidSolver::new(height, width, 0.005, 1.0 / 128.0, 0.1, bodies)
        .integration(integration::Integration::BogackiShampine)
        .linear_solver(linear_solvers::LinearSolver::ConjugateGradient {
            auxiliary: vec![0.0; width * height],
            search: vec![0.0; width * height],
            preconditioner: vec![0.0; width * height],
            a: Sparse::new(width * height),
            iterations: 600
        })
        .advection(advection::Advection::SemiLagrangian);

    // Over 10s
    for iteration in 0..500 {
        for i in 0..4 {
            print!("Step {}: ", 4 * iteration + i);
            solver.add_inflow(0.45, 0.2, 0.15, 0.03, 1.0, 0.0, 3.0);
            solver.update();
        }

        solver.to_image(1.0, &mut buffer);
        println!("Saving.");
        image::save_buffer(format!("output/{}.png", iteration), &buffer, width as u32, height as u32, image::RGBA(8)).unwrap();
    }
}