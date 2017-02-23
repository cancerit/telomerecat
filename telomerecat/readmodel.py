import pdb
import numpy as np
from  multiprocessing import Pool,freeze_support
from functools import partial

from abc import ABCMeta, abstractmethod

class Distribution(object):
    __metaclass__ = ABCMeta

    def __init__(self, params):
        self.params = params

    def __check_params__(self):
        for i in xrange(len(self.params["shape"])):
            if round(self.params["shape"][i],2) <= 0.001:
                self.params["shape"][i] = 0.001

    def copy(self):
        new_params = {"count":self.params["count"],
                      "shape":list(self.params["shape"])}

        return type(self)(new_params)

    def modify(self, param, value):
        self.params[param] = value
        self.__check_params__()

    def filter_simulated_data(self, data, lo_limit, hi_limit):
        return data
        # filtered_data = []

        # for d in data:
        #     if d < lo_limit:
        #         filtered_data.append(lo_limit)
        #     elif d > hi_limit:
        #         filtered_data.append(hi_limit)
        #     else:
        #         filtered_data.append(d)

        # return np.array(filtered_data)

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
                       "shape":[0.9,1.1,50]}
        self.distributions["atelo"] = Beta(telo_params)

        subtelo_count = observed_dist[twenty_five_percent:fifty_percent].sum()
        subtelo_params = {"count":subtelo_count,
                          "shape":[ten_percent, fifty_percent]}
        self.distributions["subtelo"] = Uniform(subtelo_params)

        nontelo_params = {"count":observed_dist[fifty_percent:].sum(),
                          "shape":[fifty_percent,ten_percent]}
        self.distributions["nontelo"] = Normal(nontelo_params)

    def copy(self):

        new_dists = {}

        for dist_name, dist in self.distributions.items():
            new_dists[dist_name] = dist.copy()

        return Solution(new_dists)

    def get_new_solutions(self, new_count, modify_dist, modify_param):

        new_solutions = []

        for _ in xrange(new_count):
            new_solution = self.copy()

            relevant_dist = new_solution.distributions[modify_dist]
            relevant_param = relevant_dist.params[modify_param]

            if modify_param == "count":
                dist_names = self.distributions.keys()

                dist_change = 0

                for dist_name, cur_dist in new_solution.distributions.items():
                    if dist_name != modify_dist:
                        cur_count = cur_dist.params["count"]

                        max_change = cur_count*.2
                        change = int(np.random.uniform(0, max_change))
                        if cur_count-change <= 0:
                            change = 0

                        dist_change += change

                        cur_dist.modify("count", cur_count - change)

                relevant_dist.modify("count", relevant_param + dist_change)

                new_solutions.append(new_solution)

            elif modify_param == "shape":
                exisiting_value = relevant_param
                new_values = [] 

                for exisiting_value in relevant_param:
                    new_values.append(np.random.normal(exisiting_value,
                                                       exisiting_value*.5))

                relevant_dist.modify(modify_param, 
                                     new_values)

                new_solutions.append(new_solution)

        return new_solutions

    def set_score(self, new_score):
        self.score = new_score

    def simulate(self, read_len):

        all_data = []
        for dist_name, params in self.distributions.items():
            sim_data = params.simulate(hi_limit=read_len)
            all_data.extend(sim_data)

        return np.array(all_data)
        
class TelomereReadModel(object):

    def __init__(self, sample_stats, observed_dist=None, read_stats=None):
        self.sample_stats = sample_stats

        if observed_dist is not None:
            self.observed_dist = observed_dist
        else:
            self.read_stats = read_stats
            mask = self.read_stats[:,0] > -1
            self.observed_data = self.read_stats[mask,0]
            self.observed_dist = self.__dist_from_data__(self.observed_data, 
                                                    self.sample_stats["read_len"])

    def __get_best_solution__(self, best_solution, solutions):
        solutions.sort(key = (lambda s: s.score))

        if solutions[0].score < best_solution.score:
            new_best_solution = solutions[0]
        else:
            new_best_solution = best_solution

        return new_best_solution

    def __score_solutions__(self, solutions):

        read_len = self.sample_stats["read_len"]

        for solution in solutions:

            simulated_data = solution.simulate(read_len)
            simulated_dist = self.__dist_from_data__(simulated_data, read_len)

            score = self.__compare_dists__(self.observed_dist, simulated_dist)
            solution.set_score(score)

    def __compare_dists__(self, observed, simulated):

        #change = ( np.log2(observed+0.0001) - np.log2(simulated+0.0001) ) **2
        dif = (observed[:40]-simulated[:40])**2
        #score = dif.sum() * change.sum()
        return int(dif.sum())

    def __dist_from_data__(self, data, read_len):
        bins = max((1, int(data.max()))) 

        dist = np.histogram(data,bins = bins)[0]
        dist_buffer = [0] * int(read_len - len(dist))
        return np.concatenate((dist, dist_buffer))

    def __print_best_solution__(self, best_solution):

        def dist_to_string(dist):
            shape_string = ",".join(["%.3f" % (s,) \
                                     for s in  dist.params["shape"]])
            return "%s|%d" % (shape_string, dist.params["count"])

        dist_string = ""
        for dist_name, dist in best_solution.distributions.items():
            dist_string += "%s: (%s) " % (dist_name[:2],dist_to_string(dist)) 

        print "%d:: %s" % (best_solution.score, dist_string)

    def to_r(self,best_solution, read_len=100):
        data = best_solution.simulate(read_len)
        dist = self.__dist_from_data__(data, read_len)
        print "sim=c(" + ",".join(["%d" % (d,) for d in dist]) + ")"

    def get_permutations(self, param_count, dist_count):
        perms = np.array(np.meshgrid(range(param_count), range(dist_count))).T
        perms = perms.reshape(-1,2).tolist()
        return perms

    def run(self):

        best_solution = Solution()
        best_solution.bootstrap(self.observed_dist)
        prev_best_solution = best_solution

        parameters = ["count","shape"]
        distributions = sorted(best_solution.distributions.keys())

        param_count = len(parameters)
        dist_count = len(distributions)

        permutations = self.get_permutations(param_count, 
                                             dist_count)

        no_change = 0
        iteration = 0
        converged = False

        while not converged:

            param_type, dist_type = permutations[iteration % \
                                                 len(permutations) ]

            modify_dist = distributions[dist_type]
            modify_param = parameters[param_type]

            new_solutions = best_solution.get_new_solutions(25, 
                                                            modify_dist, 
                                                            modify_param)

            self.__score_solutions__(new_solutions)

            prev_best_solution = best_solution
            best_solution = self.__get_best_solution__(best_solution, new_solutions) 
            #self.__print_best_solution__(best_solution)

            if best_solution.score == prev_best_solution.score:
                no_change += 1
            else:
                no_change = 0

            if no_change == 100:
                converged = True

            if iteration % 10 == 0:
                print "its:%d chg:%d scr:%d" % (iteration, 
                                                no_change, 
                                                best_solution.score)

            iteration += 1

        print "FINISHED:", best_solution.score
        return best_solution

def get_best_solution(results):
    results.sort(key = lambda s: s.score)
    return results[0]

def model_process(job, sample_stats, observed_dist):
    np.random.seed()
    model_estimator = TelomereReadModel(sample_stats, 
                                         observed_dist=observed_dist)
    best_solution = model_estimator.run()
    return best_solution

def run_model_par(sample_stats, read_stats, n_procs, N=20):

    primary_read_model = TelomereReadModel(sample_stats, 
                                           read_stats=read_stats)
    observed_dist = primary_read_model.observed_dist

    freeze_support()
    p = Pool(n_procs)
    
    model_partial = partial(model_process, 
                            sample_stats=sample_stats, 
                            observed_dist=observed_dist)
    results = p.map(model_partial,range(N))
    p.close()

    best_solution = get_best_solution(results)

    return best_solution
                

            
