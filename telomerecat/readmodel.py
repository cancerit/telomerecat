
from collections import namedtuple

class Params(object):

    def __init__(self, count, p_list):
        self.count = count
        self.p_list = plist

class Solution(object):

    def __init__(self,telo=None, 
                      subtelo=None, 
                      nontelo=None):

        self.telo = telo
        self.subtelo = subtelo
        self.nontelo = nontelo

    def bootstrap(self, observed_dist):

        read_len = len(observed_dist)
        ten_percent = int(read_len * .1)
        twenty_five_percent = int(read_len * .25)
        fifty_percent = int(read_len * .5)

        telo_count = observed_dist[:twenty_five_percent].sum()
        self.telo = Params(telo_count,[5,])
        
        subtelo_count = observed_dist[twenty_five_percent:fifty_percent].sum()
        self.subtelo = Params(subtelo_count,[twenty_five_percent,ten_percent])

        nontelo_count = observed_dist[fifty_percent:].sum()
        self.nontelo = Params(nontelo_count, [fifty_percent+ten_percent,
                                              ten_percent])

class TelomereReadModel(object):

    def __init__(self):
        pass

    def __dist_from_data__(self, data, read_len):
        dist = np.histogram(data,bins = int(data.max()) )[0]
        dist_buffer = [0] * int(read_len - len(dist))
        return np.concatenate((dist, dist_buffer))

    def model_counts(self, read_stats, sample_stats):

        converged = False

        total_reads = read_stats.shape[0]

        observed_data = read_stats[:,0]
        observed_dist = self.__dist_from_data__(observed_data, 
                                                sample_stats["read_len"])

        telo_N = observed_dist[:10].sum()
        subtelo_N = observed_dist[10:30].sum()
        nontelo_N = observed_dist[30:].sum()

        start_solution = Solution()
        start_solution.bootstrap(observed_dist)

        telo_lambda = 2

        subtelo_mu = 20
        subtelo_sigma = 5

        nontelo_mu = 60
        nontelo_sigma = 5

        while not converged:
            telo_reads = np.random.exponential(telo_lambda,size=telo_N)
            subtelo_reads = np.random.normal(subtelo_mu,subtelo_sigma,subtelo_N)
            nontelo_reads = np.random.normal(nontelo_mu,nontelo_sigma,nontelo_N)

            sim_data = np.concatenate((telo_reads, 
                                       subtelo_reads, 
                                       nontelo_reads))

            sim_dist = self.__dist_from_data__(sim_data, 
                                               sample_stats["read_len"])
            pdb.set_trace()