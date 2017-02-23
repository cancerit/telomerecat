import pdb
import numpy as np
from collections import namedtuple

from abc import ABCMeta, abstractmethod

class Distribution(object):
    __metaclass__ = ABCMeta

    def __init__(self, params):
        self.params = params

    def copy(self):
        new_params = {"count":self.params["count"],
                      "shape":list(self.params["shape"])}

        return type(self)(new_params)

    def modify(self, param, index, value):
        if param == "count":
            self.params[param] = value
        else:
            self.params[param][index] = value

    def filter_simulated_data(self, data, lo_limit, hi_limit):
        filtered_data = []

        for d in data:
            if d < lo_limit:
                filtered_data.append(lo_limit)
            elif d > hi_limit:
                filtered_data.append(hi_limit)
            else:
                filtered_data.append(d)

        return np.array(filtered_data)

    @abstractmethod
    def simulate(self):
        pass

    @abstractmethod
    def bootstrap(self):
        pass

class Exponential(Distribution):

    def __init__(self,params):
        super(Exponential, self).__init__(params)

    def simulate(self, lo_limit = 0, hi_limit = 100):
        data = np.random.exponential(self.params["shape"][0],
                                     self.params["count"])
        return self.filter_simulated_data(data, lo_limit, hi_limit)

    def bootstrap(self):
        param1 = np.random.uniform(0,5)
        self.params["shape"] = [param1]

class Beta(Distribution):

    def __init__(self,params):
        super(Beta, self).__init__(params)

    def simulate(self, lo_limit = 0, hi_limit = 100):
        data = np.random.beta(self.params["shape"][0],
                              self.params["shape"][1],
                              self.params["count"])
        data = data * self.params["shape"][2]

        return self.filter_simulated_data(data, lo_limit, hi_limit)

    def bootstrap(self):
        param1 = np.random.uniform(0,8)
        param2 = np.random.uniform(0,8)
        param3 = np.random.uniform(0,100)
        self.params["shape"] = [param1, param2, param3]

class Uniform(Distribution):

    def __init__(self,params):
        super(Uniform, self).__init__(params)

    def simulate(self, lo_limit = 0, hi_limit = 100):
        data = np.random.uniform(self.params["shape"][0],
                                 self.params["shape"][1],
                                 self.params["count"])
        return self.filter_simulated_data(data, lo_limit, hi_limit) 

    def bootstrap(self):
        param1 = np.random.uniform(0,99)
        param2 = np.random.uniform(param1,100)
        self.params["shape"] = [param1, param2]

class Normal(Distribution):

    def __init__(self,params):
        super(Normal, self).__init__(params)

    def simulate(self, lo_limit = 0, hi_limit = 100):
        data = np.random.normal(self.params["shape"][0],
                                self.params["shape"][1],
                                self.params["count"])

        return self.filter_simulated_data(data, lo_limit, hi_limit)

    def bootstrap(self):
        param1 = np.random.uniform(0,100)
        param2 = np.random.uniform(1,10)
        self.params["shape"] = [param1, param2]

class Solution(object):

    def __init__(self,distributions = None):

        if distributions is None:
            self.distributions = {}
        else:
            self.distributions = distributions
        self.score = float("Inf")

    def bootstrap(self, observed_dist):

        read_len = len(observed_dist)
        ten_percent = int(read_len * .1)
        twenty_five_percent = int(read_len * .25)
        fifty_percent = int(read_len * .5)

        telo_params = {"count":observed_dist[:twenty_five_percent].sum(),
                       "shape":[3,]}
        self.distributions["telo"] = Exponential(telo_params)

        subtelo_count = observed_dist[twenty_five_percent:fifty_percent].sum()
        subtelo_params = {"count":subtelo_count,
                          "shape":[twenty_five_percent, ten_percent]}
        self.distributions["subtelo"] = Normal(subtelo_params)

        #self.subtelo.bootstrap()
        # subtelo_count = observed_dist[twenty_five_percent:fifty_percent].sum()
        # subtelo_params = {"count":subtelo_count,
        #                   "shape":[ten_percent, fifty_percent]}
        # self.subtelo = Uniform(subtelo_params)

        # nontelo_params = {"count":observed_dist[fifty_percent:].sum(),
        #                   "shape":[fifty_percent+ten_percent, ten_percent]}
        # self.nontelo = Normal(nontelo_params)

        # nontelo_params = {"count":observed_dist[fifty_percent:].sum(),
        #                   "shape":[5,2,70]}
        # self.nontelo = Beta(nontelo_params)
        # self.nontelo.bootstrap()

    def copy(self):

        new_telo = self.telo.copy()
        new_subtelo = self.subtelo.copy()
        new_nontelo = self.nontelo.copy()

        return Solution(new_telo,
                        new_subtelo,
                        new_nontelo)

    def get_new_solutions(self, new_count, modify_dist, modify_param):

        new_solutions = []

        for _ in xrange(new_count):
            new_solution = self.copy()

            relevant_dist = getattr(new_solution, modify_dist)
            relevant_param = relevant_dist.params[modify_param]

            if modify_param == "count":
                dist_names = self.distributions.keys()

                dist_change = 0

                for dist_name in dist_names:
                    if dist_name != modify_dist:
                        cur_dist = getattr(new_solution, dist_name)
                        cur_count = cur_dist.params["count"]

                        max_change = cur_count*.2
                        change = int(np.random.uniform(0, max_change))
                        if cur_count-change <= 0:
                            change = 0

                        dist_change += change

                        cur_dist.modify("count",-1, cur_count - change)

                relevant_dist.modify("count", -1, relevant_param + dist_change)

                new_solutions.append(new_solution)

            elif modify_param == "shape":
                modify_index = np.random.choice(xrange(len(relevant_param)))

                exisiting_value = relevant_param[modify_index]
                new_value = (np.random.normal(exisiting_value,
                                              exisiting_value*.2))

                relevant_dist.modify(modify_param, 
                                     modify_index, 
                                     new_value)

                new_solutions.append(new_solution)

        return new_solutions

    def set_score(self, new_score):

        self.score = new_score

    def simulate(self, read_len):
        telo_data = self.telo.simulate(hi_limit=read_len)
        subtelo_data = self.subtelo.simulate(hi_limit=read_len)
        nontelo_data = self.nontelo.simulate(hi_limit=read_len)

        return np.concatenate((telo_data, subtelo_data, nontelo_data))
        
class TelomereReadModel(object):

    def __init__(self, read_stats, sample_stats):
        self.read_stats = read_stats
        self.sample_stats = sample_stats

    def __get_best_solution__(self, best_solution, solutions):
        solutions.sort(key = (lambda s: s.score))

        if solutions[0].score < best_solution.score:
            new_best_solution = solutions[0]
        else:
            new_best_solution = best_solution

        return new_best_solution

    def __score_solutions__(self, solutions, observed_dist):

        read_len = self.sample_stats["read_len"]

        for solution in solutions:

            simulated_data = solution.simulate(read_len)
            simulated_dist = self.__dist_from_data__(simulated_data, read_len)

            score = self.__compare_dists__(observed_dist, simulated_dist)
            solution.set_score(score)

    def __compare_dists__(self, observed, simulated):

        dif = ( np.log2(observed+0.0001) - np.log2(simulated+0.0001) ) **2
        score = dif.sum()
        return score

    def __dist_from_data__(self, data, read_len):
        try:
            bins = max((1, int(data.max()))) 
        except ValueError:
            pdb.set_trace()

        dist = np.histogram(data,bins = bins)[0]
        dist_buffer = [0] * int(read_len - len(dist))
        return np.concatenate((dist, dist_buffer))

    def __print_best_solution__(self, best_solution):

        def dist_to_string(dist):
            shape_string = ",".join(["%.3f" % (s,) \
                                     for s in  dist.params["shape"]])
            return "%s|%d" % (shape_string, dist.params["count"])

        print "%d:: t:E(%s) s:N(%s) n:N(%s)" % \
                    (best_solution.score,
                     dist_to_string(best_solution.telo),
                     dist_to_string(best_solution.subtelo),
                     dist_to_string(best_solution.nontelo))

    def to_r(self,best_solution, read_len=100):
        data = best_solution.simulate(read_len)
        dist = self.__dist_from_data__(data, read_len)
        print "sim=c(" + ",".join(["%d" % (d,) for d in dist]) + ")"

    def get_permutations(self, param_count, dist_count):
        perms = np.array(np.meshgrid(range(param_count), range(dist_count))).T
        perms = perms.reshape(-1,2).tolist()
        return perms.tolist()

    def run(self):

        converged = False

        mask = self.read_stats[:,0]<40
        observed_data = self.read_stats[mask,0]
        observed_dist = self.__dist_from_data__(observed_data, 
                                                self.sample_stats["read_len"])

        best_solution = Solution()
        best_solution.bootstrap(observed_dist)

        iteration = 0

        parameters = ["count","shape"]
        distributions = best_solution.distributions.keys()

        param_count = len(parameters)
        dist_count = len(distributions)

        permutations = self.get_permutations(param_count, 
                                             dist_count)

        while not converged:

            param_type, dist_type = permutations[iteration % \
                                                 (param_count+dist_count) ]

            modify_dist = distributions[dist_type]
            modify_param = parameters[param_type]

            new_solutions = best_solution.get_new_solutions(20, 
                                                            modify_dist, 
                                                            modify_param)

            self.__score_solutions__(new_solutions, observed_dist)

            best_solution = self.__get_best_solution__(best_solution, new_solutions) 

            #[self.__print_best_solution__(s) for s in new_solutions]

            self.__print_best_solution__(best_solution)

            iteration += 1

            if iteration % 100 == 0:
                self.to_r(best_solution)
                pdb.set_trace()
            #print iteration % 2, iteration % 3
            



            






