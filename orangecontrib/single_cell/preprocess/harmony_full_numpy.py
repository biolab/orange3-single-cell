import numpy as np

def _l2norm_cols(X, eps=1e-12):
    nrm = np.linalg.norm(X, axis=0, keepdims=True)
    nrm = np.maximum(nrm, eps)
    return X / nrm

def _one_hot_labels(labels):
    labels = np.asarray(labels)
    levels, inv = np.unique(labels, return_inverse=True)
    B = len(levels)
    N = labels.shape[0]
    Phi = np.zeros((B, N), dtype=float)
    Phi[inv, np.arange(N)] = 1.0
    counts = Phi.sum(1)
    return Phi, levels, counts

def _one_hot_multi(vars_list):
    if not isinstance(vars_list, (list, tuple)):
        vars_list = [vars_list]
    Phis, levels_all, counts_all = [], [], []
    for v in vars_list:
        Phi, lev, cnt = _one_hot_labels(v)
        Phis.append(Phi)
        levels_all.append(lev)
        counts_all.append(cnt)
    Phi = np.concatenate(Phis, axis=0)
    counts_concat = np.concatenate(counts_all)
    return Phi, levels_all, counts_concat

def _expand_theta(theta, counts_per_var):
    counts_per_var = np.asarray(counts_per_var)
    if np.isscalar(theta):
        return np.repeat(float(theta), counts_per_var.sum())
    theta = np.asarray(theta, dtype=float)
    if theta.ndim == 1 and theta.size == counts_per_var.size:
        return np.concatenate([np.repeat(theta[j], counts_per_var[j]) for j in range(len(counts_per_var))])
    if theta.ndim == 1 and theta.size == counts_per_var.sum():
        return theta
    raise ValueError("theta: unknown shape.")

def _kpp_cosine_init(Z_cos_T, K, rng):
    N, d = Z_cos_T.shape
    centers = np.empty((K, d), dtype=float)
    idx = rng.integers(0, N)
    centers[0] = Z_cos_T[idx]
    D = 1.0 - (Z_cos_T @ centers[0][:, None]).squeeze()
    D = np.maximum(D, 0.0)
    for k in range(1, K):
        probs = D / (D.sum() + 1e-12)
        idx = rng.choice(N, p=probs)
        centers[k] = Z_cos_T[idx]
        dots = Z_cos_T @ centers[k][:, None]
        D = np.minimum(D, 1.0 - dots.squeeze())
        D = np.maximum(D, 0.0)
    centers /= np.linalg.norm(centers, axis=1, keepdims=True).clip(min=1e-12)
    return centers

def _compute_dist_cos(Y_T, Z_cos):
    return 2.0 * (1.0 - (Y_T @ Z_cos))

def _update_R_blocks(R, dist_mat, sigma, Phi, Pr_b, theta_level, block_size, rng, E, O):
    K, N = R.shape
    B = Phi.shape[0]

    score = -dist_mat / sigma[:, None]
    score -= np.max(score, axis=0, keepdims=True)
    base = np.exp(score)

    order = np.arange(N)
    rng.shuffle(order)
    n_blocks = int(np.ceil(1.0 / block_size))
    blocks = np.array_split(order, n_blocks)

    for b in blocks:
        if b.size == 0:
            continue

        mass = R[:, b].sum(axis=1)
        E -= np.outer(mass, Pr_b)
        O -= R[:, b] @ Phi[:, b].T

        penalty_kb = np.power((E + 1.0) / (O + 1.0), theta_level)
        penalty_cells = penalty_kb @ Phi[:, b]
        R[:, b] = base[:, b] * penalty_cells
        colsum = np.sum(R[:, b], axis=0, keepdims=True).clip(min=1e-12)
        R[:, b] /= colsum

        mass_new = R[:, b].sum(axis=1)
        E += np.outer(mass_new, Pr_b)
        O += R[:, b] @ Phi[:, b].T

    return R, E, O

def _moe_correct_ridge_numpy(Z_orig, Z_corr, Z_cos, R, Phi_moe, lamb):
    d, N = Z_orig.shape
    K = R.shape[0]
    Z_corr[:] = Z_orig
    for k in range(K):
        Phi_Rk = Phi_moe * R[k, :][None, :]
        X = Phi_Rk @ Phi_moe.T
        X += lamb
        RHS = Phi_Rk @ Z_orig.T
        W = np.linalg.solve(X, RHS)
        W[0, :] = 0.0
        Z_corr -= (W.T @ Phi_Rk)
    Z_cos[:] = _l2norm_cols(Z_corr, 1e-12)
    return Z_corr, Z_cos

def _safe_entropy(R):
    out = R * np.log(R, where=(R>0), out=np.zeros_like(R))
    return out

def _compute_objective_full(R, dist, sigma, theta_level, E, O, Phi):
    kmeans_error = float(np.sum(R * dist))
    ent = float(np.sum(_safe_entropy(R) * sigma[:, None]))
    y = np.tile(theta_level[None, :], (R.shape[0], 1))
    z = np.log((O + 1.0) / (E + 1.0))
    w = (y * z) @ Phi
    cross = float(np.sum((R * sigma[:, None]) * w))
    total = kmeans_error + ent + cross
    return total, kmeans_error, ent, cross

def harmony_numpy_fun(
    X_raw,               
    batches,
    n_pcs=50,             
    K=None,
    sigma=0.1,
    theta=1.0,
    lamb=1.0,
    tau=0.0,
    block_size=0.05,
    max_iter_harmony=10,
    max_iter_kmeans=20,
    eps_kmeans=1e-5,
    eps_harmony=1e-4,
    random_state=0,
    print_every=1
):
    rng = np.random.default_rng(random_state)

    X = np.asarray(X_raw, dtype=float)
    Xc = X - X.mean(axis=0, keepdims=True)          # center
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    Z = (U[:, :n_pcs] * S[:n_pcs])                  # PCA embedding

    N, d = Z.shape

    Phi, levels_all, counts_concat = _one_hot_multi(batches)
    B = Phi.shape[0]
    Pr_b = Phi.sum(axis=1) / N

    counts_per_var = np.array([len(lev) for lev in levels_all])
    theta_level = _expand_theta(theta, counts_per_var)

    if tau > 0 and K is not None:
        theta_level = theta_level * (1.0 - np.exp(-(counts_concat / (K * tau))**2))

    if np.isscalar(lamb):
        lamb_vec = np.full(B + 1, float(lamb))
        lamb_vec[0] = 0.0
        lamb_mat = np.diag(lamb_vec)
    else:
        lamb_vec = np.asarray(lamb, dtype=float)
        lamb_mat = np.diag(lamb_vec)

    Phi_moe = np.vstack([np.ones((1, N)), Phi])

    Z_T = Z.T.copy()
    colmax = np.max(np.abs(Z_T), axis=0, keepdims=True).clip(min=1e-12)
    Z_cos = _l2norm_cols(Z_T / colmax, 1e-12)
    Z_corr = Z_T.copy()

    if K is None:
        K = int(min(round(N / 30.0), 100))

    if np.isscalar(sigma):
        sigma = np.full(K, float(sigma))
    else:
        sigma = np.asarray(sigma, dtype=float)

    Y = _kpp_cosine_init(Z_cos.T, K, rng)
    Y /= np.linalg.norm(Y, axis=1, keepdims=True).clip(min=1e-12)

    dist = _compute_dist_cos(Y, Z_cos)
    score = -dist / sigma[:, None]
    score -= np.max(score, axis=0, keepdims=True)
    R = np.exp(score)
    R /= np.sum(R, axis=0, keepdims=True)

    E = np.outer(R.sum(axis=1), Pr_b)
    O = R @ Phi.T

    obj_hist = []
    total, km, ent, cross = _compute_objective_full(R, dist, sigma, theta_level, E, O, Phi)
    obj_hist.append(total)
    print(f"[init] Obj={total:.6e} | kmeans={km:.6e} | entropy={ent:.6e} | cross={cross:.6e}")

    for it in range(1, max_iter_harmony + 1):
        for it_k in range(1, max_iter_kmeans + 1):

            Y = (R @ Z_cos.T)
            Y /= np.linalg.norm(Y, axis=1, keepdims=True).clip(min=1e-12)

            dist = _compute_dist_cos(Y, Z_cos)
            R, E, O = _update_R_blocks(R, dist, sigma, Phi, Pr_b, theta_level,
                                       block_size, rng, E, O)

        Z_corr, Z_cos = _moe_correct_ridge_numpy(Z_T, Z_corr, Z_cos, R, Phi_moe, lamb_mat)

        total, km, ent, cross = _compute_objective_full(R, dist, sigma, theta_level, E, O, Phi)
        obj_hist.append(total)

        if it % print_every == 0:
            print(f"[outer {it:02d}] Obj={total:.6e} | kmeans={km:.6e} | entropy={ent:.6e} | cross={cross:.6e}")

        if len(obj_hist) >= 2:
            rel = (obj_hist[-2] - obj_hist[-1]) / (abs(obj_hist[-2]) + 1e-12)
            if rel < eps_harmony:
                print(f"[stop] Convergence achieved: Δrel={rel:.3e} < eps_harmony={eps_harmony}")
                break

    return Z_corr.T   
