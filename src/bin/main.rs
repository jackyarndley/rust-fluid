use rust_fluid::fluid_solver::FluidSolver;
use rust_fluid::{integration, linear_solvers, advection};
use rust_fluid::util::sparse::Sparse;

extern crate image;

fn main() {
    let width = 256;
    let height = 128;

    let mut buffer = vec![0u8; width * height * 4];

    let mut solver = FluidSolver::new(height, width, 0.005, 1.0 / 128.0, 1.0)
        .integration(integration::Integration::BogackiShampine)
        .linear_solver(linear_solvers::LinearSolver::ConjugateGradient {
            auxiliary: vec![0.0; width * height],
            search: vec![0.0; width * height],
            preconditioner: vec![0.0; width * height],
            a: Sparse::new(width * height),
            iterations: 600
        })
        .advection(advection::Advection::SemiLagrangian);

    solver.u_velocity.add_inflow(20, 20, 40, 40, 5.0);
    solver.v_velocity.add_inflow(20, 20, 40, 40, 0.0);
    solver.density.add_inflow(20, 20, 40, 40, 2.0);

    for iteration in 0..500 {
        for i in 0..4 {
            print!("Iteration {}. ", 4 * iteration + i);
            solver.update();

            solver.u_velocity.add_inflow(20, 20, 40, 40, 5.0);
            solver.v_velocity.add_inflow(20, 20, 40, 40, 0.0);
            solver.density.add_inflow(20, 20, 40, 40, 2.0);
        }

        solver.to_image(2.0, &mut buffer);
        println!("Saving image...");
        image::save_buffer(format!("output/{}.png", iteration), &buffer, width as u32, height as u32, image::RGBA(8)).unwrap();
    }
}