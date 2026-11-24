use diffsol::{
    Bdf, BdfState, FaerLU, FaerMat, FaerVec, OdeBuilder, OdeSolverMethod, OdeSolverState, Vector,
};

#[test]
fn native_bdf_relative_error_scale_uses_the_new_zero_crossing_state() {
    let h: f64 = 0.001;
    let rtol: f64 = 1e-6;
    let atol: f64 = 1e-10;
    let problem = OdeBuilder::<FaerMat<f64>>::new()
        .t0(0.0)
        .h0(h)
        .rtol(rtol)
        .atol([atol])
        .rhs_implicit(
            |_x, _p, t, out| out[0] = -1000.0 + 0.001 * t,
            |_x, _p, _t, _seed, out| out[0] = 0.0,
        )
        .init(|_p, _t, out| out[0] = 1.0, 1)
        .build()
        .unwrap();
    let state = BdfState::<FaerVec<f64>>::new_without_initialise(&problem).unwrap();
    let mut method: Bdf<'_, _, _, FaerMat<f64>> =
        Bdf::new(&problem, state, FaerLU::<f64>::default()).unwrap();
    method.step().unwrap();
    let t = method.state().t;
    let actual = method.state().y.get_index(0);
    let exact = 1.0 - 1000.0 * t + 0.0005 * t * t;
    let budget = atol + rtol * exact.abs();
    assert!(
        (actual - exact).abs() <= budget,
        "accepted y={actual} at{t}, exact={exact}, error={} exceeds new-state budget={budget}",
        (actual - exact).abs()
    );
    assert!(
        method.get_statistics().number_of_error_test_failures >= 1,
        "the original trial must be rejected at its new zero-crossing coordinate"
    );
    assert!(t < h);
}
