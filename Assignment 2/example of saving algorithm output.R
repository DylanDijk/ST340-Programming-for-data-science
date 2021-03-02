load('data/20newsgroups.rdata') 
source('data/em_mixture_bernoullis.R')

load("poo.rds")
poo = em_mix_bernoulli(xs[1:4,], K  = 4)       
save(poo, file = "poo.rds")



xs = documents
out_Q1bi = em_mix_bernoulli(xs, K  = 4)  
save(out_Q1bi, file = "output_of_EM_mix_bernoulli.rds")


load("output_of_EM_mix_bernoulli.rds")

xs = documents
out_mu_4_6 = em_mix_bernoulli_mu_4_6(xs, K  = 4)
save(out_Q1bi, out_mu_4_6,file = "output_of_EM_mix_bernoulli.rds")

xs = documents
out_W_prop = em_mix_bernoulli_W_prop(xs, K  = 4)
save(out_W_prop, out_Q1bi, out_mu_4_6,file = "output_of_EM_mix_bernoulli.rds")


xs = documents
out_W_mu_prop  =  em_mix_bernoulli_W_mu_prop(xs, K = 4)
save(out_W_mu_prop, out_W_prop, out_Q1bi, out_mu_4_6,file = "output_of_EM_mix_bernoulli.rds")


xs = documents[sample(nrow(documents), size = 1000), ]

em_mix_bernoulli_test(xs, K = 4)
