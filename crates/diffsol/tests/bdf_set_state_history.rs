use diffsol::{
    Bdf, BdfState, FaerLU, FaerMat, FaerVec, OdeBuilder, OdeSolverMethod, OdeSolverState,
};

#[test]
fn replaced_bdf_state_starts_a_new_equal_step_history() {
    let problem = OdeBuilder::<FaerMat<f64>>::new()
        .t0(0.0)
        .h0(0.001)
        .rtol(1e-6)
        .atol([1e-10])
        .rhs_implicit(
            |_x, _p, _t, out| out[0] = 0.0,
            |_x, _p, _t, _seed, out| out[0] = 0.0,
        )
        .init(|_p, _t, out| out[0] = 1.0, 1)
        .build()
        .unwrap();
    let state = BdfState::<FaerVec<f64>>::new_without_initialise(&problem).unwrap();
    let mut method: Bdf<'_, _, _, FaerMat<f64>> =
        Bdf::new(&problem, state, FaerLU::<f64>::default()).unwrap();
    method.step().unwrap();
    assert_eq!(method.order(), 1);
    let time = method.state().t;
    let mut fresh = BdfState::<FaerVec<f64>>::new_without_initialise(&problem).unwrap();
    *fresh.as_mut().t = time;
    fresh.set_problem(&problem).unwrap();
    method.set_state(fresh);
    method.step().unwrap();
    assert_eq!(
        method.order(),
        1,
        "first step of replaced history must not reuse equal steps from the old history"
    );
}
