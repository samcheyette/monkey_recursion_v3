
import random
import copy
#from context import *
from cleanStevesData_byParticipant2 import *
import math
from scipy.stats import beta

matching_set = [('(', ')'), ('[', ']')]
open_set = ["(", "["]
closed_set = [")", "]"]
all_p = ['[', '(', ')', ']']

types = ['O', 'C', 'M'] + all_p



def get_hyps(prims, length=4):
	if length == 1:
		return prims

	else:
		new_hs = []
		for h in get_hyps(prims, length-1):
			for p in prims:
				new_hs.append(h + p)
		return copy.deepcopy(new_hs)

def get_hyp_gen(hyp, available,
				 open_available, 
				 closed_available, 
				 sofar=""):


	def find_match(m, available):
	    for r in available:
	        if (m, r) in matching_set or (r, m) in matching_set:
	            return r
	    return None
	
	if len(hyp) == 0:
		return {"":1.0}

	else:
		h = hyp[0]
		poss = {}
		if h in "([])" and h in available:
			poss[h] = 1.0

		elif h == "O" and len(open_available)>0:
			for o in open_available:
				poss[o] = 1/float(len(open_available))

		elif h == "C" and len(closed_available)>0:
			for o in closed_available:
				poss[o] = 1/float(len(closed_available))

		elif h == "M" and len(closed_available) > 0:
			for s in sofar[::-1]:
				match = find_match(s, available)
				if (match in closed_available):
					poss[match] = 1.0
					break

		if len(poss.keys()) == 0:
			poss["*"] = 1.0


		ret_dcts = {}
		for key in poss:
			new_av = copy.copy(available)
			new_opav = copy.copy(open_available)
			new_clav = copy.copy(closed_available)
			new_sofar = sofar + key
			if key in new_av:
				new_av.remove(key)
				if key in new_opav:
					new_opav.remove(key)
				else:
					new_clav.remove(key)

			from_here = get_hyp_gen(hyp[1:], 
									new_av,
									new_opav,
									new_clav,
									new_sofar)
			prob_key = poss[key]
			for k in from_here:
				prob_here =from_here[k]
				if key + k not in ret_dcts:
					ret_dcts[key + k] = 0.0
				ret_dcts[key+k] += prob_here * prob_key

		return ret_dcts


def get_hyps_gen(hyps):
	r_dct = {}
	for hyp in hyps:
		hyp_gen = get_hyp_gen(hyp, 
			copy.copy(all_p),
			copy.copy(open_set), 
			copy.copy(closed_set))

		r_dct[hyp] = copy.deepcopy(hyp_gen)
	return r_dct


def filter_hyps(out_hyps, thresh=2.0):
	keep = {}
	for h in out_hyps:
		badness = 0.0
		for k in out_hyps[h]:
			badness += float(k.count("*")) * out_hyps[h][k]



		if badness < thresh:
			keep[h] = out_hyps[h]
	return keep

def remove_dups(out_hyps):

	seen = set()
	new_dct = dict()
	for i in xrange(4, -1, -1):
		for o in sorted(out_hyps.keys())[::-1]:
			if ((o.count("(") +
				 o.count("[") + 
				 o.count("]") + 
				 o.count(")")) >= i):

				#not necessary and may not always work...
				#just for clean look.......
				o = o[:2] + o[2:].replace("CM", "CC")
				o = o[:2] + o[2:].replace("MC", "MM")

				lst_form =[]
				for x in out_hyps[o]:
					lst_form.append((x, out_hyps[o][x]))
				lst_form = tuple(sorted(lst_form))
				if lst_form not in seen:
					seen.add(lst_form)
					new_dct[o] = out_hyps[o]
	return new_dct


def prob_d_given_h(hyp_gen, data, sm=1e-10):
	n = sum([data[d] for d in data])

	p_h = 0.0
	hyp_gen_use = copy.deepcopy(hyp_gen)
	sum_val = float(sum([data[k] for k in data]))
	z = 0.0
	for d in data:
		if d not in hyp_gen_use:
			hyp_gen_use[d] = 0.0
		hyp_gen_use[d] += sm/sum_val

	for h in hyp_gen_use:
		z += hyp_gen_use[h]

	for d in data:
		#p_h += math.log(sm)
		p_h += (math.log(hyp_gen_use[d]) - math.log(z)) * (data[d]) 
		#p_h += math.log(hyp_gen[d]) * data[d] 

	return p_h

def normalize(lst, z=None):
	z = sum(lst)
	assert(z > 0.0)
	for l in lst:
		new_lst.append(l/float(z))
	return z

def flatten_theta(theta, out_hyps, okeys):

	new_out = {}
	for t in xrange(len(theta)):

		for k in out_hyps[okeys[t]]:
			if k not in new_out:
				new_out[k] = 0.0
			new_out[k] += out_hyps[okeys[t]][k] * theta[t]

	return copy.deepcopy(new_out)



def output_full_model(data, out_file, hyp_distr):
	o = "who, which, P, alpha, CE, CR, TL, OT, corr\n"

	ce = ["([])", "[()]"]
	cr = ["([)]", "[(])"]
	tr = ["[]()", "()[]"]

	for d in data:
		who = d[0]
		alpha = d[1]
		beta = d[2]
		for b in beta:
			o += "%s, %s, %f, %f, " % (who,
										 b, 
										 beta[b], 
										 alpha)


			p_ce = 0.0
			p_cr = 0.0
			p_tr = 0.0
			p_ot = 0.0

			for h in hyp_distr[b]:
				if h in ce:
					p_ce += hyp_distr[b][h]
				elif h in cr:
					p_cr += hyp_distr[b][h]
				elif h in tr:
					p_tr += hyp_distr[b][h]
				else:
					p_ot += hyp_distr[b][h]

			corr=False
			if p_ce == 1:
				corr=True

			o += "%f, %f, %f, %f, %s\n" % (p_ce, p_cr, p_tr, p_ot, str(corr))



	wr = open(out_file, "w+")
	wr.write(o)
	wr.close()