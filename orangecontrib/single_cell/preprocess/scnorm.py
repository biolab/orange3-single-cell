import numpy as np
import scipy as sp
import itertools as it
from Orange.regression import Learner, Model, SklLearner
from statsmodels.regression.quantile_regression import QuantReg
from sklearn.preprocessing import PolynomialFeatures


class QuantRegLearner(Learner):
    """Generate a quantile regression modeland learn a prediction model

    Parameters
    ----------
    degree : float in [0, 1]
        target quantile
    preprocessors : List[Preprocessor]
        preprocessors to be applied on the data before learning
    """
    name = 'quantile regression'
    preprocessors = SklLearner.preprocessors

    def __init__(self, quantile=0.5, preprocessors=None):
        super().__init__(preprocessors=preprocessors)
        self.quantile = quantile

    def fit(self, X, Y, W):
        qr = QuantReg(Y, X)
        model = qr.fit(self.quantile)
        params = model.params
        return QuantRegModel(model, params)


class QuantRegModel(Model):
    def __init__(self, model, params):
        self.model = model
        self.params = params

    def predict(self, X):
        return self.model.predict(X)

    def __str__(self):
        return 'QuantRegModel {}'.format(self.model)


# TODO: Split genes into K groups using K-medioids and merge groups smaller than 100 genes

def scnorm(data, subgroup_frac = 0.25, K = 10):
    """
    ScNorm normalization method.

    :param data: Orange data Table with cells in rows and genes in colums.
    :param subgroup_frac: Subgroup fraction.
    :param K: number of groups.
    :return:
    """

    # Copy data to be transformed
    new_data = data.copy()

    # Fixed algorithm parameters
    q_range = np.linspace(0.05, 0.95, 19)
    degree_range = np.arange(1, 7)

    # Calculate sequencing depth for each cell
    n, m = data.X.shape
    D = data.X.sum(axis=1).reshape((n, 1))

    # Fit a quantile regression model for each gene based on sequencing depth
    model_slopes = dict()
    median_expressions = dict()
    for var in data.domain.attributes:
        Y = data.get_column_view(var)[0].reshape((n, 1))
        rows = np.nonzero(Y)[0]
        if len(rows) > 1:
            Dp = PolynomialFeatures(degree=1).fit_transform(np.log(D[rows]))
            depth_model = QuantRegLearner().fit(Dp, np.log(Y[rows]), W=None)
            model_slopes[var] = depth_model.params[1]
            median_expressions[var] = np.median(Y[rows])


    # Partitioning of genes into equally-sized groups
    valid_vars = list(model_slopes.keys())
    exps = np.array(list(map(lambda v: median_expressions[v], valid_vars)))
    percentiles = np.array(list(map(lambda e: sp.stats.percentileofscore(exps, e), exps)))
    groups = np.array(K * percentiles / 100.0, dtype=int)
    groups[groups == K] = K - 1

    # Non-zero slopes
    slopes = np.array(list(map(lambda v: model_slopes[v], valid_vars)))

    # Process each group independently
    for k in range(K):
        print("Processing group %d ..." % k)
        group_slopes = slopes[groups == k]
        group_median = np.median(group_slopes)
        group_inxs = np.argsort(np.absolute(group_slopes - group_median))[:max(2, int(subgroup_frac * len(group_slopes)))]
        group_vars = [valid_vars[i] for i in np.where(groups == k)[0][group_inxs]]
        group_counts = data[:, group_vars].X

        # Fill submodel data
        nz = np.nonzero(group_counts)
        group_D = np.zeros((len(nz[0]), 1))
        group_Y = np.zeros((len(nz[0]), 1))
        for i, (r, c) in enumerate(zip(*nz)):
            group_D[i] = D[r]
            group_Y[i] = group_counts[r, c]

        # Fit submodels for different quantiles and select best model
        submodels = dict()
        q_best = degree_best = q_distance = float("inf")
        for q, degree in it.product(q_range, degree_range):
            Dp = PolynomialFeatures(degree=degree).fit_transform(np.log(group_D))
            submodels[q, degree] = QuantRegLearner(quantile=q).fit(Dp,
                                                                   np.log(group_Y),
                                                                   W=None)
            curr_dist = np.absolute(group_median - submodels[q, degree].params[1])
            if curr_dist < q_distance:
                q_distance = curr_dist
                q_best = q
                degree_best = degree
        print("Best model parameters (q=%.2f, d=%d)" % (q_best, degree_best))

        # Select best model and use normalized values
        Dp = PolynomialFeatures(degree=degree_best).fit_transform(np.log(group_D))
        log_predicted_expression = submodels[q_best, degree_best].predict(Dp)
        gene_expression = sp.stats.scoreatpercentile(group_Y, 100 * q_best)
        size_factors = np.exp(log_predicted_expression) / np.exp(gene_expression)
        new_data[:, group_vars].X[nz] = new_data[:, group_vars].X[nz] / size_factors

    # Return data with a normalized matrix
    return new_data



if __name__ == "__main__":
    from Orange.data.table import Table
    import matplotlib.pyplot as plt
    test_file = "/Users/martin/Dev/data/singlecell/zhang2017/tables/matrix_counts_sample.tab"
    data = Table(test_file)
    cols = np.where((data.X > 0).sum(axis=0) > 1)[0]
    data = data[:, cols]
    data.X *= np.random.rand(*data.X.shape)

    # Normalize original data
    new_data = scnorm(data, K=1)
    x = np.log(data.X.ravel() + 1)
    y = np.log(new_data.X.ravel() + 1)

    plt.figure()
    plt.plot(x, y, ".")
    plt.show()