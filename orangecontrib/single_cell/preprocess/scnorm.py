import numpy as np
import scipy as sp
from Orange.regression import Learner, Model, SklLearner
from statsmodels.regression.quantile_regression import QuantReg


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



def scnorm(data):
    """
    ScNorm normalization method
    :param data: Orange data Table with cells in rows and genes in colums.
    :return:
    """

    # Calculate sequencing depth for each cell
    K = 10
    n, m = data.X.shape
    D = data.X.sum(axis=1).reshape((n, 1))

    # Fit a quantile regression model for each gene based on sequencing depth
    depth_models = dict()
    model_slopes = dict()
    median_expressions = dict()
    for var in data.domain.attributes:
        Y = data.get_column_view(var)[0].reshape((n, 1))
        rows = np.nonzero(Y)[0]
        if len(rows) > 1:
            depth_models[var] = QuantRegLearner().fit(np.log(D[rows]),
                                                      np.log(Y[rows]),
                                                      W=None)
            model_slopes[var] = depth_models[var].params[0]
            median_expressions[var] = np.median(Y[rows])


    # Partitioning of genes into equally-sized groups
    # TODO: Split genes into K groups using K-medioids
    G = np.zeros((m, K))
    exps = np.array(list(map(lambda v: median_expressions[v], data.domain.attributes)))
    slopes = np.array(list(map(lambda v: slopes[v], data.domain.attributes)))
    percentiles = np.array(list(map(lambda e: sp.stats.percentileofscore(exps, e), exps)))
    groups = np.array(K * percentiles / 100.0, dtype=int)
    groups[groups == K] = K - 1
    for i, j in enumerate(groups): G[i, j] = 1

    # Debug: plot
    si = np.argmax(slopes)
    var = data.domain.attributes[si]
    Y = data.get_column_view(var)[0].reshape((n, 1))
    rows = np.nonzero(Y)[0]
    logy = depth_models[var].predict(np.log(D[rows]))
    plt.plot(np.log(D[rows]), np.log(Y[rows]), "k.")
    plt.plot(np.log(D[rows]), logy, "b-")
    plt.plot(np.log(D[rows]), slopes[si] * np.log(D[rows]), "--", color="gray")
    plt.xlabel("Log sequencing depth")
    plt.ylabel("Log expression")
    plt.title(var.name)
    plt.show()

    # TODO: Process each group to compute its size factors
    Sf = np.zeros((n, K))  # Size factors to be determined

    # Return a normalized matrix
    Xnew = data.X * data.Sf.dot(G.T)
    return Table.from_numpy(X=Xnew,
                            Y=data.Y,
                            metas=data.metas,
                            W=data.W,
                            domain=data.domain)



if __name__ == "__main__":
    from Orange.data.table import Table
    import matplotlib.pyplot as plt
    test_file = "/Users/martin/Dev/data/singlecell/zhang2017/tables/matrix_counts_sample.tab"
    data = Table(test_file)
    cols = np.where((data.X > 0).sum(axis=0) > 1)[0]
    data = data[:, cols]
    data.X *= np.random.rand(*data.X.shape)

    learner = QuantRegLearner()
    X, y = data.X[:, [1,3]], data.X[:, 5:6]
    model = learner.fit(X, y, W=None)