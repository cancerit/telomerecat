import time
import pdb
import numpy as np
from  multiprocessing import Pool,freeze_support
from functools import partial

from abc import ABCMeta, abstractmethod


class Distribution(object):
    __metaclass__ = ABCMeta

    def __init__(self, params, observed_dist, lo, hi):
        if params is None:
            self.params = {}
            self.params["count"] = observed_dist[lo:hi].sum()
            self.bootstrap(observed_dist, lo, hi)
        else:
            self.params = params

    def __check_params__(self):
        for i in xrange(len(self.params["shape"])):
            if round(self.params["shape"][i], 2) <= 0.001:
                self.params["shape"][i] = 0.001

    def copy(self):
        new_params = {"count": self.params["count"],
                      "shape": list(self.params["shape"])}

        return type(self)(params=new_params)

    def modify(self, param, value):
        self.params[param] = value
        self.__check_params__()

    def dist_to_dat(self, dist, offset=0):
        dat = []
        for i, x in enumerate(dist):
            dat.extend([i + offset, ] * x)
        return np.array(dat)

    @abstractmethod
    def simulate(self):
        pass

    @abstractmethod
    def bootstrap(self, observed_dist, lo, hi):
        pass


class Exponential(Distribution):

    def __init__(self, params=None, observed_dist=None, lo=None, hi=None):
        super(Exponential, self).__init__(params, observed_dist, lo, hi)

    def simulate(self, lo_limit=0, hi_limit=100):
        data = np.random.exponential(self.params["shape"][0],
                                     self.params["count"])
        return data

    def bootstrap(self, observed_dist, lo, hi):
        param1 = np.random.uniform(0, 5)
        self.params["shape"] = [param1]


class Beta(Distribution):

    def __init__(self, params=None, observed_dist=None, lo=None, hi=None):
        super(Beta, self).__init__(params, observed_dist, lo, hi)

    def simulate(self, lo_limit=0, hi_limit=100):
        data = np.random.beta(self.params["shape"][0],
                              self.params["shape"][1],
                              self.params["count"])
        data = data * self.params["shape"][2]

        return data

    def bootstrap(self, observed_dist, lo, hi):
        self.params["shape"] = self.get_start_params(observed_dist, lo, hi)

    def get_start_params(self, observed_dist, lo, hi):
        def equation(x_bar, v_bar):
            return ((x_bar * (1 - x_bar)) / v_bar) - 1

        dat = self.dist_to_dat(observed_dist[lo:hi], offset=lo)
        max_obs = dat.max()
        dat = dat / (max_obs + .001)

        dat_mean = dat.mean()
        dat_var = dat.var()

        alpha = dat_mean * equation(dat_mean, dat_var)
        beta = (1 - dat_mean) * equation(dat_mean, dat_var)

        return [alpha, beta, max_obs]


class Uniform(Distribution):

    def __init__(self, params=None, observed_dist=None, lo=None, hi=None):
        super(Uniform, self).__init__(params, observed_dist, lo, hi)

    def simulate(self, lo_limit=0, hi_limit=100):
        data = np.random.uniform(self.params["shape"][0],
                                 self.params["shape"][1],
                                 self.params["count"])
        return data

    def bootstrap(self, observed_dist, lo, hi):
        self.params["shape"] = [lo, hi]


class Normal(Distribution):

    def __init__(self, params=None, observed_dist=None, lo=None, hi=None):
        super(Normal, self).__init__(params, observed_dist, lo, hi)

    def simulate(self, lo_limit=0, hi_limit=100):
        data = np.random.normal(self.params["shape"][0],
                                self.params["shape"][1],
                                self.params["count"])

        return data

    def bootstrap(self, observed_dist, lo, hi):
        dat = self.dist_to_dat(observed_dist[lo:hi], offset=lo)
        self.params["shape"] = [dat.mean(), dat.std()]


class Solution(object):

    def __init__(self, read_len, distributions=None):

        self.read_len = read_len

        if distributions is None:
            self.distributions = {}
        else:
            self.distributions = distributions
        self.score = float("Inf")

    def bootstrap(self, observed_dist):

        read_len = len(observed_dist)
        ten_percent = int(read_len * .1)
        twenty_percent = int(read_len * .2)
        fifty_percent = int(read_len * .5)

        self.distributions["atelo"] = Beta(observed_dist=observed_dist,
                                           lo=0,
                                           hi=twenty_percent)

        self.distributions["subtelo"] = Uniform(observed_dist=observed_dist,
                                                lo=twenty_percent,
                                                hi=fifty_percent)

        self.distributions["nontelo"] = Normal(observed_dist=observed_dist,
                                               lo=fifty_percent,
                                               hi=fifty_percent + ten_percent)

    def copy(self):

        new_dists = {}

        for dist_name, dist in self.distributions.items():
            new_dists[dist_name] = dist.copy()

        return Solution(self.read_len, new_dists)

    def new_count_params(self, modify_dist):
        modify_param = "count"

        relevant_dist = self.distributions[modify_dist]
        relevant_param = relevant_dist.params[modify_param]

        dist_change = 0
        for dist_name, cur_dist in self.distributions.items():
            if dist_name != modify_dist:
                cur_count = cur_dist.params[modify_param]

                max_change = cur_count * .2
                change = int(np.random.uniform(0, max_change))
                if cur_count - change <= 0:
                    change = 0

                dist_change += change

                cur_dist.modify("count", cur_count - change)

        relevant_dist.modify("count", relevant_param + dist_change)

    def new_shape_params(self, modify_dist):
        modify_param = "shape"

        relevant_dist = self.distributions[modify_dist]
        relevant_param = relevant_dist.params[modify_param]

        exisiting_value = relevant_param
        new_values = []

        for exisiting_value in relevant_param:
            new_values.append(np.random.normal(exisiting_value,
                                               exisiting_value))
        relevant_dist.modify(modify_param,
                             new_values)

    def get_new_solutions(self, new_count, modify_dist, modify_param):

        new_solutions = []

        for _ in xrange(new_count):
            new_solution = self.copy()
            if modify_param == "count":
                new_solution.new_count_params(modify_dist)
            elif modify_param == "shape":
                new_solution.new_shape_params(modify_dist)
            new_solutions.append(new_solution)

        return new_solutions

    def set_score(self, new_score):
        self.score = new_score

    def simulate(self, exclude_dists=[]):

        all_data = []
        for dist_name, params in self.distributions.items():
            if dist_name not in exclude_dists:
                sim_data = params.simulate(hi_limit=self.read_len)
                all_data.extend(sim_data)

        return np.array(all_data)

    def get_dist(self, exclude_dists=[]):
        data = self.simulate(exclude_dists)
        return self.data_to_dist(data)

    def data_to_dist(self, data):
        bins = max((1, int(data.max())))
        dist = np.histogram(data, bins=bins)[0]
        dist_buffer = [0] * int(self.read_len - len(dist))
        return np.concatenate((dist, dist_buffer))

    def model_to_r(self):
        dist = self.get_dist()
        self.print_r(dist)

    def print_r(self, dist):
        print "sim=c(" + ",".join(["%d" % (d,) for d in dist]) + ")"

    def print_params(self):
        def dist_to_string(dist):
            shape_string = ",".join(["%.3f" % (s,)
                                     for s in dist.params["shape"]])
            return "%s|%d" % (shape_string, dist.params["count"])

        dist_string = ""
        for dist_name, dist in self.distributions.items():
            dist_string += "%s: (%s) " % (dist_name[:2], dist_to_string(dist))

        print "%d:: %s" % (self.score, dist_string)

class TelomereReadModel(object):

    def __init__(self, sample_stats,
                       observed_dist=None,
                       read_stats=None,
                       job_id=0):

        self.sample_stats = sample_stats
        self.job_id = job_id

        if observed_dist is not None:
            self.observed_dist = observed_dist
        else:
            self.read_stats = read_stats
            mask = self.read_stats[:, 0] > -1
            self.observed_data = self.read_stats[mask, 0]

            self.observed_dist = \
                    Solution(self.sample_stats["read_len"]).\
                            data_to_dist(self.observed_data)

    def __get_best_solution__(self, best_solution, solutions):
        solutions.sort(key=(lambda s: s.score))

        if solutions[0].score < best_solution.score:
            new_best_solution = solutions[0]
        else:
            new_best_solution = best_solution

        return new_best_solution

    def __score_solutions__(self, solutions):

        for solution in solutions:
            simulated_dist = solution.get_dist()
            score = self.__compare_dists__(self.observed_dist, simulated_dist)
            solution.set_score(score)

    def __compare_dists__(self, observed, simulated):

        # change = ( np.log2(observed+0.0001) - np.log2(simulated+0.0001) ) **2
        dif = (observed[:40] - simulated[:40])**2
        # score = dif.sum() * change.sum()
        return int(dif.sum())

    def get_permutations(self, param_count, dist_count):
        perms = np.array(np.meshgrid(range(param_count), range(dist_count))).T
        perms = perms.reshape(-1, 2).tolist()
        return perms

    def run(self):

        best_solution = Solution(self.sample_stats["read_len"])
        best_solution.bootstrap(self.observed_dist)
        prev_best_solution = best_solution

        parameters = ["count", "shape"]
        distributions = sorted(best_solution.distributions.keys())

        param_count = len(parameters)
        dist_count = len(distributions)

        permutations = self.get_permutations(param_count,
                                             dist_count)

        no_change = 0
        iteration = 0
        converged = False

        start = time.time()

        while not converged:

            param_type, dist_type = permutations[iteration % len(permutations)]

            modify_dist = distributions[dist_type]
            modify_param = parameters[param_type]

            new_solutions = best_solution.get_new_solutions(10,
                                                            modify_dist,
                                                            modify_param)

            self.__score_solutions__(new_solutions)

            prev_best_solution = best_solution
            best_solution = self.__get_best_solution__(best_solution,
                                                       new_solutions)
            
            if best_solution.score == prev_best_solution.score:
                no_change += 1
            else:
                no_change = 0

            if no_change == 100:
                converged = True

            # if iteration % 10 == 0:
            #     best_solution.print_params()
            #     print "%d:: its:%d chg:%d scr:%d" % (self.job_id,
            #                                          iteration,
            #                                          no_change,
            #                                          best_solution.score)

            iteration += 1

        print "JOB %d:: SCORE: %d TIME:%d" % \
                  (self.job_id,
                   best_solution.score,
                   time.time() - start,)

        return best_solution


def get_best_solution(results):
    results.sort(key=lambda s: s.score)
    return results[0]


def model_process(job, sample_stats, observed_dist):
    np.random.seed()
    model_estimator = TelomereReadModel(sample_stats,
                                        observed_dist=observed_dist,
                                        job_id=job)
    best_solution = model_estimator.run()
    return best_solution


def get_read_model(sample_stats, read_stats, n_procs, N=20):

    primary_read_model = TelomereReadModel(sample_stats,
                                           read_stats=read_stats)
    observed_dist = primary_read_model.observed_dist

    freeze_support()
    p = Pool(n_procs)

    model_partial = partial(model_process,
                            sample_stats=sample_stats,
                            observed_dist=observed_dist)
    results = p.map(model_partial, range(N))
    p.close()

    best_solution = get_best_solution(results)

    return best_solution



