# 04/21/2021

# an updated version of inference utilities
# for simulation & inference purposes

library(tidyverse)
library(foreach)
library(doParallel)
library(Matrix)
registerDoParallel()

# 0. to prepare stuff
# (1) a function to generate samples from Exp(rate) truncated at (a,b)
sam_texp <- function(n, rate, a=0, b=Inf){
  quans = runif(n)
  sams = -(log(exp(-rate*a)-quans*(exp(-rate*a)-exp(-rate*b))))/rate
  return(sams)
}

# (2.1) function to get adjmat up to time point t
get_adjmat <- function(G0, events, t){
  adjmat = G0
  net_events = events %>% 
    filter(event %in% c(3:8), time < t)
  
  if(nrow(net_events) > 0){
    for(r in 1:nrow(net_events)){
      z = net_events$event[r]
      p1 = net_events$per1[r]
      p2 = net_events$per2[r]
      if(z %in% c(3:5)){
        adjmat[p1,p2] = 1
        adjmat[p2,p1] = 1
      }else{
        adjmat[p1,p2] = 0
        adjmat[p2,p1] = 0
      }
    }
  }
  
  adjmat
}

# (2.2) a function to save adjmat structure at some report time points
get_adjmat_report <- function(G0, events, times){
  
  adjmat = G0
  net_events = events %>% filter(event %in% c(3:8))
  
  res = list()
  
  if(times[1]==0){
    lb = times[1]
    res[[as.character(lb)]] = adjmat
  }
  
  if(nrow(net_events)==0){
    for(i in c(1:length(times))){
      lb = times[i]
      res[[as.character(lb)]] = adjmat
    }
  }else{
    t = 0; r = 1
    for(i in c(1:length(times))){
      lb = times[i]
      while((t<lb) & (r <= nrow(net_events))){
        t = net_events$time[r]
        z = net_events$event[r]
        p1 = net_events$per1[r]
        p2 = net_events$per2[r]
        
        if(t >= lb){
          res[[as.character(lb)]] = adjmat
        }
        
        if(z %in% c(3:5)){
          adjmat[p1,p2] = 1
          adjmat[p2,p1] = 1
        }else{
          adjmat[p1,p2] = 0
          adjmat[p2,p1] = 0
        }
        r = r+1
      }
      
      if(t < lb){
        res[[as.character(lb)]] = adjmat
      }
    }
  }
  
  res
}

# # try this out
# times = seq(0,56,by=7)
# G_report = get_adjmat_report(G0, events, times)
# # test out to make sure things are correct!
# for(t in times){
#   adj_t = get_adjmat(G0, events, t)
#   cat(all(adj_t == G_report[[as.character(t)]]),"\n")
#   cat(which(adj_t != G_report[[as.character(t)]], arr.ind = TRUE),'\n')
# }

# events %>% filter(event %in% c(3:8), 
#                   (per1==161 & per2==145)|(per1==145 & per2==161))


# (3) a function to generate missing exposure &/ recovery times
#     and produce status reports

miss_data <- function(events, G0, I0, interval = 7, 
                      miss_recov_prop = 1, miss_expo_prop = 1){
  # `interval`: the regular period of status report
  # `miss_recov_prop`: the proportion/probability of missing recovery times
  # `miss_expo_prop`: the proportion/probability of missing exposure times
  
  N = nrow(G0)
  ## 10/19/2020: epid[I0] should be -1 (E)
  epid = rep(0,N); epid[I0] = -1
  tmax = max(events$time)
  report.times = seq(from = 0, to = (tmax %/% interval + 1) * interval, by = interval)
  report = epid
  
  # select all epi events
  # events.epi = events[events$event %in% c(1,2),]
  
  recov.ind = NULL
  expo.ind = NULL
  for(ix in c(2:length(report.times))){
    lb = report.times[ix-1]; ub = report.times[ix]
    st = min(which(events$time > lb)); en = max(which(events$time <= ub))
    if(st == Inf){ break }
    for(r in st:en){
      z = events$event[r]
      p = events$per1[r]
      if(z==1){
        # exposure
        epid[p] = -1
        expo.ind = c(expo.ind, r)
      }else if(z==2){
        # recovery
        epid[p] = -2
        recov.ind = c(recov.ind, r)
      }else if(z %in% c(9,10)){
        # manifestation
        epid[p] = ifelse(z==9, 1, 2)
      }
    }
    # status report
    report = rbind(report, epid)
  }
  
  # records to remove in recovery times and exposure times
  remove_recov = as.logical(rbinom(length(recov.ind), 1, miss_recov_prop))
  remove_expo = as.logical(rbinom(length(expo.ind), 1, miss_expo_prop))
  to_remove = c(recov.ind[remove_recov], expo.ind[remove_expo])
  events = events[-to_remove,]
  
  # return new dataset and status report
  row.names(report) = as.character(report.times)
  
  # if("truth" %in% names(dat)){
  #   return(list(G0=G0, I0=I0, report.times=report.times, 
  #               events=events, report = report, 
  #               truth = dat$truth))
  # }else{
  #   cat("True parameters unavailable for this dataset!\n")
  #   return(list(G0=G0, I0=I0, report.times=report.times, 
  #               events=events, report = report))
  # }
  
  list(G0=G0, I0=I0, report.times=report.times, 
       events=events, report = report)
}

# # try this
# miss_dat = miss_data(events, G0, I0)


# (4) a function to read in simulated data (generated in Python)
load_data <- function(root, dirname){
  # root: has to be the full path
  # dirname: the folder name
  
  #fpath = paste0(root,dirname)
  fpath = file.path(root,dirname)
  
  # events
  events = read_csv(file.path(fpath,'dat.csv'), col_types = cols())
  names(events) = c('time','event','per1','per2')
  ## need to re-label people with 1-indexing (rather than 0-indexing)
  events$per1 = events$per1 + 1
  events$per2 = events$per2 + 1
  
  # G0 and X
  G0 = as.matrix(read_delim(file.path(fpath,'G0.txt'), delim=" ", 
                            col_names = FALSE, col_types = cols()))
  colnames(G0) = NULL
  ## use parse matrix to save some space
  G0 = as(G0, "dgTMatrix")
  X = as.matrix(read_delim(file.path(fpath,'X.txt'), delim=" ", 
                           col_names = FALSE, col_types = cols()))
  colnames(X) = NULL
  
  # the params, I0 and change_times
  params = read_csv(file.path(fpath,'params.csv'), col_types = cols())
  I0 = params$I0[1] + 1 # need to add 1 too
  stage_change = params$stage_change
  
  # also record exp(eta)
  truth = list(beta = params$beta[1], eta = params$eta[1], exp_eta = exp(params$eta[1]),
               gamma = params$gamma[1], phi = params$phi[1], p_s = params$p_s[1])
  
  # the network params
  net_params = read_csv(file.path(fpath,'net_params.csv'), col_types = cols())
  truth$alpha = net_params$alpha
  truth$omega = net_params$omega
  
  # the regression coefficients
  regr_params = read_csv(file.path(fpath,'regr_params.csv'), col_types = cols())
  truth$b_S = regr_params$b_S
  truth$b_alpha = regr_params$b_alpha
  truth$b_omega = regr_params$b_omega
  
  # put everything together
  dat = list(events = events, G0 = G0, I0 = I0, stage_change = stage_change, 
             X = X, truth = truth)
  
  dat
}

## try this out
# dat = load_data('~/Documents/Research_and_References/Hetero_EpiNet_2020/', 'ex1')


#####
# I. for MLEs

# (1) summarize epi table for each person i
# function to summarize epi data for person i
summarize_epi <- function(i, G0_i, I0, events){
  
  # init epid vector to record disease status
  N = length(G0_i)
  epid = rep(0,N)
  epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G0_i
  # make sure (i,i) entry is 0
  nei[i] = 0
  
  # data storage for 
  # individual epidemic info
  epi_tab = numeric(6)
  names(epi_tab) = c("Ia_expo", "Is_expo", "Ia_ti", "Is_ti",
                     "sick_time", "latent_time")
  
  # get tmax
  tmax = max(events$time)
  
  # epidemic events for i (E, Ia/Is, R times)
  events_i = events %>% filter(per1 == i, event %in% c(1,9,10,2))
  # calculate sick time and latent time if i was ever infected
  if(nrow(events_i)>0){
    #expo_time = events_i$time[1]
    
    # get latent time
    ill_ind = which(events_i$event %in% c(9,10))
    if(length(ill_ind) > 0){
      epi_tab["latent_time"] = events_i$time[ill_ind] - events_i$time[1]
      # then get ill time
      recov_ind = which(events_i$event==2)
      if(length(recov_ind) > 0){
        epi_tab['sick_time'] = events_i$time[recov_ind] - events_i$time[ill_ind]
      }else{
        epi_tab['sick_time'] = tmax - events_i$time[ill_ind]
      }
    }else{
      epi_tab["latent_time"] = tmax - events_i$time[1]
    }
    
    # set "tmax" to exposed time if ever infected...
    tmax = events_i$time[1]
  }
  
  # get all events that are related to i (both epi and net)
  # AND all the manisfestation and recovery events (event %in% c(2,9,10))
  # ONLY need events that are before Exposure time (if no exposure then tmax)
  events = events %>% filter(time <= tmax) %>% 
    filter(per1 == i | per2 == i | event %in% c(2,9,10))
  
  # if there is nothing, then NO exposure for i (nobody ever manisfested)
  if(nrow(events) > 0){
    # go through the events one by one
    t_pre = 0
    for(r in 1:nrow(events)){
      z = events$event[r]
      t_cur = events$time[r]
      
      # obtain cumulated time of exposure up to now
      epi_tab['Ia_expo'] = epi_tab['Ia_expo'] + sum(nei * (epid==1)) * (t_cur - t_pre)
      epi_tab['Is_expo'] = epi_tab['Is_expo'] + sum(nei * (epid==2)) * (t_cur - t_pre)
      
      # then process changes to the system
      if (z==1){
        # exposure
        p1 = events$per1[r]
        epid[p1] = -1
        
        # get local neighborhood at time of exposure
        if(p1 == i){
          epi_tab['Ia_ti'] = sum(nei * (epid==1))
          epi_tab['Is_ti'] = sum(nei * (epid==2))
          
          # some info
          # taken out!
          #cat('when ',p1,' got infected, had ', epi_tab['Ia_ti'],
          #    'Ia contacts and ', epi_tab['Is_ti'], 'Is contacts.\n')
        }
      }else if (z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        p1 = events$per1[r]
        epid[p1] = ifelse(z==9, 1, 2)
      }else if (z==2){
        # recovery
        p1 = events$per1[r]
        epid[p1] = -2 # changed coding -2=R
      }else{
        # some edge stuff
        p1 = events$per1[r]
        p2 = events$per2[r]
        
        if(z %in% c(3:5)){
          # reconnection
          if(p1 == i){
            nei[p2] = 1
          }else if(p2 == i){
            nei[p1] = 1
          }
        }else if(z %in% c(6:8)){
          # disconnection
          if(p1 == i){
            nei[p2] = 0
          }else if(p2 == i){
            nei[p1] = 0
          }
        }
      }
      # set t_pre to t_cur
      t_pre = t_cur
    }
  }
  
  epi_tab
}


# (2) summarize connection/disconnection history for (i,j) pair

# a helper function to get the "link type"
# used to index the c_ij or d_ij vector for i,j network events
# given sub-vector of epid[i,j] and "stage" (binary: 0 or 1)
get_sub_type <- function(epid){
  res = numeric(3)
  
  total = sum(epid)
  if(total == 0){
    # both H
    res[1] = 1
  }else if(total == 1){
    # H and I
    res[2] = 1
  }else{
    # both I
    res[3] = 1
  }
  
  res
}

# function to process info for i,j pairs
# UPDATE 10/18/2020: also keep track of C and D counts (to sum up later on)
# output length: 7 + 7 + 6 + 6
summarize_ij <- function(i,j, G0_ij, I0, events, stage_change){
  # tmax
  tmax = max(events$time)
  
  # select all events related to i,j
  events = events %>% filter((per1 == i & per2==j) | (per1 == j & per2 == i) | 
                               (per1==i & is.na(per2)) | (per1==j & is.na(per2)))
  
  # get st and en of NPI period
  st = stage_change[1]; en = stage_change[2]
  
  # storage for C and D counts
  c_counts = rep(0,6)
  d_counts = rep(0,6)
  type_phase = c('HH0','HI0','II0','HH1','HI1','II1')
  names(c_counts) = type_phase
  names(d_counts) = type_phase
  
  # if there is no event related to i or j...
  # i,j have stayed H
  # and G_ij == G0_ij throughout
  if(nrow(events)==0){
    if(G0_ij==1){
      # they've stayed connected
      c_ij = rep(0,7)
      d_ij = rep(0,7)
      d_ij[2] = (st-0) + (tmax-en)
      d_ij[5] = en - st
    }else{
      # they've stayed disconnected
      d_ij = rep(0,7)
      c_ij = rep(0,7)
      c_ij[2] = (st-0) + (tmax-en)
      c_ij[5] = en - st
    }
  }else{
    # else: something happened for them
    epid = numeric(2)
    G_ij = G0_ij
    
    # # separate two stages
    # events_T0 = events %>% filter(time <= st | time >= en)
    # events_T1 = events %>% filter(time > st & time < en)
    
    c_ij = rep(0,7)
    d_ij = rep(0,7)
    
    t_pre = 0
    
    # !
    # attach a final fake row 
    # to account for the dwell time till the end
    fake_row  = c(tmax, 666, -1, -1)
    events = rbind(events, fake_row)
    
    # go through events one by one
    for(r in 1:nrow(events)){
      
      # output something
      #print(c_ij)
      #print(d_ij)
      
      z = events$event[r]
      t_cur = events$time[r]
      
      link_type = get_sub_type(epid)
      # accumulate risks (dwell times)
      if(t_cur < st | t_pre > en){
        # operate under T0
        if(G_ij == 0){
          # spent disconnected
          c_ij[2:4] = c_ij[2:4] + link_type * (t_cur - t_pre)
        }else{
          # spent connected
          d_ij[2:4] = d_ij[2:4] + link_type * (t_cur - t_pre)
        }
      }else if(t_pre > st & t_cur < en){
        # operate under T1
        if(G_ij == 0){
          # spent disconnected
          c_ij[5:7] = c_ij[5:7] + link_type * (t_cur - t_pre)
        }else{
          # spent connected
          d_ij[5:7] = d_ij[5:7] + link_type * (t_cur - t_pre)
        }
      }else if(t_pre < st & t_cur > en){
        # covers the T1 period: some T0 + entire T1 + some T0
        if(G_ij == 0){
          # spent disconnected
          c_ij[2:4] = c_ij[2:4] + link_type * (st - t_pre + t_cur - en)
          c_ij[5:7] = c_ij[5:7] + link_type * (en - st)
        }else{
          # spent connected
          d_ij[2:4] = d_ij[2:4] + link_type * (st - t_pre + t_cur - en)
          d_ij[5:7] = d_ij[5:7] + link_type * (en - st)
        }
      }else{
        # straddling T0 and T1
        # both T0 and T1 need update
        
        ## get change point first:
        cp = ifelse(t_cur < en, st, en)
        
        #cat('straddle case: change point =', cp, 't_pre = ', t_pre, 't_cur = ', t_cur, '\n')
        
        if(G_ij == 0){
          # spent disconnected
          if(cp==st){
            c_ij[2:4] = c_ij[2:4] + link_type * (cp - t_pre)
            c_ij[5:7] = c_ij[5:7] + link_type * (t_cur - cp)
          }else{
            c_ij[5:7] = c_ij[5:7] + link_type * (cp - t_pre)
            c_ij[2:4] = c_ij[2:4] + link_type * (t_cur - cp)
          }
        }else{
          # spent connected
          if(cp==st){
            d_ij[2:4] = d_ij[2:4] + link_type * (cp - t_pre)
            d_ij[5:7] = d_ij[5:7] + link_type * (t_cur - cp)
          }else{
            d_ij[5:7] = d_ij[5:7] + link_type * (cp - t_pre)
            d_ij[2:4] = d_ij[2:4] + link_type * (t_cur - cp)
          }
        }
      }
      
      if(z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        per1 = events$per1[r]
        if(per1==i){
          epid[1] = 1
        }else if(per1==j){
          epid[2] = 1
        }
      }else if(z == 2){
        # recovery: becoming H again
        per1 = events$per1[r]
        if(per1==i){
          epid[1] = 0
        }else if(per1==j){
          epid[2] = 0
        }
      }else if(z %in% c(3:5)){
        # connection
        G_ij = 1
        c_ij[1] = c_ij[1] + 1
        # also update c_count
        if(t_cur < st | t_cur > en){
          # operate under T0
          c_counts[1:3] = c_counts[1:3] + link_type
        }else{
          # operate under T1
          c_counts[4:6] = c_counts[4:6] + link_type
        }
      }else if(z %in% c(6:8)){
        # disconnection
        G_ij = 0
        d_ij[1] = d_ij[1] + 1
        # also update d_count
        if(t_cur < st | t_cur > en){
          # operate under T0
          d_counts[1:3] = d_counts[1:3] + link_type
        }else{
          # operate under T1
          d_counts[4:6] = d_counts[4:6] + link_type
        }
      }
      
      t_pre = t_cur
    }
    
    
    # # finally, deal with the status through the end
    # if(t_cur < st){
    #   # last event before start of T1: some T0 and then T1 (st to en) and then T0
    # }else if(t_cur > en){
    #   # last event after end of T1: some T0 (t_cur)
    # }else{
    #   # last event during T1: some T1 (t_cur to en) and then T0
    # }
    
  }
  
  # return as a longgg vector
  # UPDATE 10/18/2020: extend this with the c_count and d_count vectors as well
  # length: 7 + 7 + 6 + 6
  c(c_ij, d_ij, c_counts, d_counts)
}

# # try it out
# i=39,j=43
# i=37;j=39
# G0_ij = G0[i,j]
# summ_ij = summarize_ij(i, j, G0_ij, I0, events, c(5,30))
# summ_ij
# 
# check it against "ground truth"
# summ$net_c_table[get_vec_index(c(i,j),N),]
# summ$net_d_table[get_vec_index(c(i,j),N),]

## a helper function to generate i, j index sequences
get_ij_seq <- function(N){
  num.E = N * (N-1)/2
  I = numeric(num.E)
  J = numeric(num.E)
  
  counter = 0
  for(i in c(1:(N-1))){
    # append: (i, i+1), (i,i+2), ..., (i, N)
    re = N - i
    I[(counter+1):(counter+re)] = rep(i, re)
    J[(counter+1):(counter+re)] = c((i+1):N)
    
    counter = counter + re
  }
  
  list(I=I, J=J)
}


# (3) everything together: function to obtain summary statistics/tables
summarize_events2 <- function(G0, I0, events, stage_change){
  
  # get basics
  N = nrow(G0)
  
  # get event counts
  nE = sum(events$event==1)
  nR = sum(events$event==2)
  nIa = sum(events$event==9)
  nIs = sum(events$event==10)
  #nI = nIa + nIs
  
  # get st and en of NPI period
  st = stage_change[1]; en = stage_change[2]
  
  # separate two stages
  events_T0 = events %>% filter(time <= st | time >= en)
  events_T1 = events %>% filter(time > st & time < en)
  
  # obtain epi table
  epi_tables = foreach(i=1:N, .combine = 'rbind') %dopar% {
    G0_i = G0[i,]
    summarize_epi(i, G0_i, I0, events)
  }
  
  # obtain net tables
  IJ = get_ij_seq(N)
  I = IJ$I; J = IJ$J
  net_tables = foreach(i=I, j=J, .combine = 'rbind') %dopar% {
    G0_ij = G0[i,j]
    summarize_ij(i, j, G0_ij, I0, events, stage_change)
  }
  
  # get C counts and D counts
  type_phase = c('HH0','HI0','II0','HH1','HI1','II1')
  #C_all = numeric(6); D_all = numeric(6)
  C_all = colSums(net_tables[,15:20])
  D_all = colSums(net_tables[,21:26])
  # for(type in c(3:5)){
  #   # connections before change
  #   C_all[type-2] = sum(events_T0$event==type)
  #   # connections after change
  #   C_all[type+1] = sum(events_T1$event==type)
  # }
  # for(type in c(6:8)){
  #   # disconnections before change
  #   D_all[type-5] = sum(events_T0$event==type)
  #   # disconnections
  #   D_all[type-2] = sum(events_T1$event==type)
  # }
  names(C_all) = type_phase
  names(D_all) = type_phase
  
  # get net_c_table and net_d_table
  net_tables = net_tables[,1:14]
  net_c_table = net_tables[,1:7]
  net_d_table = net_tables[,8:14]
  
  # data frame column names
  epi_table = as.data.frame(epi_tables)
  names(epi_table) = c('Ia_expo','Is_expo','Ia_ti','Is_ti',
                       'sick_time','latent_time')
  net_c_table = as.data.frame(net_c_table)
  names(net_c_table) = c('Nc',type_phase)
  net_d_table = as.data.frame(net_d_table)
  names(net_d_table) = c('Nd',type_phase)
  
  # return summary results
  list(counts = c(nE=nE, nIa=nIa, nIs=nIs, nR=nR),
       C_all = C_all, D_all = D_all,
       epi_table = epi_table, 
       net_c_table = net_c_table,
       net_d_table = net_d_table, I0=I0)
}

## try it out
## and time it
# st = Sys.time()
# summ200 = summarize_events2(G0, I0, events, change_times)
# cat(Sys.time() - st)
# ~ 10 seconds on my laptop - not great, but not bad either


# (4): function to iteratively solve for MLEs

## (4.1) a function that sums up X_i and X_j given p=(i,j)
sum_covariates <- function(p, X){
  # make sure that i < j
  if(p[1] > p[2]){
    p = p[2:1]
  }
  i = p[1]; j = p[2]
  
  X[i,]+X[j,]
}

solve_MLE <- function(summaries, X, maxIter=10, tol=1e-4, initEta = 1, 
                      hack0=FALSE, hackb_S = FALSE,
                      hack_eta = FALSE, hack_phi = FALSE,
                      true_Eta = NULL, true_phi = NULL, true_b_S = c(0,1)){
  
  # 10/31/2020:
  # add a hack0 mode: 
  # adjust all b_* to rep(0,p) to focus on main parameters
  
  # 11/06/2020:
  # add hack_eta_phi mode
  # fix eta and phi to truth and only focus on beta!
  # also added true_Eta and true_phi as parameters if need to hack
  
  # 12/01/2020
  # add hackb_S mode
  # set b_S to the truth and see if we can get beta and eta right
  # ALSO try to restrict exp_eta > 1
  
  # 01/26/2021
  # modify to "hack_eta" and "hack_phi" mode
  
  # X has to be a dataframe with column names!!
  
  # wrangle X into needed shape of Xij (for network part)
  # (original X: n by p design matrix)
  # (want XX: n*(n-1)/2 by p design matrix)
  N = nrow(X)
  pairs = which(upper.tri(matrix(nrow=N, ncol=N)), arr.ind = T)
  XX = apply(pairs, 1, sum_covariates, X)
  XX = t(XX)
  
  
  # extract summary statistics
  nI = summaries$counts['nE'] # take nE = nI
  nIs = summaries$counts['nIs']
  nIa = summaries$counts['nIa']
  nR = summaries$counts['nR']
  C_all = summaries$C_all
  D_all = summaries$D_all
  epi_table = summaries$epi_table
  net_c_table = summaries$net_c_table
  net_d_table = summaries$net_d_table
  I0 = summaries$I0
  
  # 0. get manifestation rate phi and symp. prob. p_s
  phi = (nIs + nIa)/sum(epi_table$latent_time)
  p_s = nIs/(nIs + nIa)
  
  if(hack_phi){
    phi = true_phi
  }
  
  # 1. get recovery rate gamma
  gamma = nR/sum(epi_table$sick_time)
  
  # 2. get parameters for the infection side: beta, eta, b_S
  eta = initEta
  beta = nI/(sum(epi_table$Ia_expo) + eta * sum(epi_table$Is_expo))
  
  # 10/19/20: change a little bit to make sure things can run...
  # get rid of I0: no exposure on I0 at all
  #got_infected = which(epi_table$latent_time > 0)
  
  # 10/24/20 fix: only get those who had infec neighbors at infec time!
  got_infected = which((epi_table$Ia_ti + epi_table$Is_ti) > 0 & epi_table$latent_time > 0)
  got_infected = got_infected[got_infected != I0]
  Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
  Y_infec[I0] = 0
  
  for(it in 1:maxIter){
    # 2.1 estimate b_S
    Expo = (epi_table$Ia_expo + eta * epi_table$Is_expo) * beta
    has_expo = which(Expo > 0) # check I0 not here!
    ## set Y_infec to the original thing before re-subsetting
    Y_infec = (epi_table$latent_time > 0) %>% as.numeric()
    Y_infec[I0] = 0
    Y_infec = Y_infec[has_expo]
    X_infec = X[has_expo,]
    Expo = Expo[has_expo]
    infec_poi = glm(Y_infec~X_infec-1+offset(log(Expo)), family = poisson())
    b_S = infec_poi$coefficients %>% as.numeric()
    
    # hack it to zero if...
    if(hack0){
      b_S = rep(0, length(b_S))
    }
    
    # (12/01/2020)
    # hack b_S to true value if ...
    if(hackb_S){
      b_S = true_b_S
    }
    
    # 2.2 estimate beta
    beta_old = beta
    indiv_effects = exp(c(X_infec %*% b_S))
    beta = nI/sum(indiv_effects * Expo/beta_old)
    
    # 2.3 estimate eta
    llfunc = function(eta){
      beta * sum(indiv_effects * (epi_table$Ia_expo + eta * epi_table$Is_expo)[has_expo]) -
        sum(log((epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected]))
    }
    ## 10/18/2020: de-bugged - previously this was not correct
    llgrr = function(eta){
      beta * sum(indiv_effects * epi_table$Is_expo[has_expo]) -
        sum(epi_table$Is_ti[got_infected]/(epi_table$Ia_ti + epi_table$Is_ti * eta)[got_infected])
    }
    eta_old = eta
    
    ## output some info
    #cat('beta=',beta, 'eta=', eta, '\nfn(eta)=', llfunc(eta),
    #   'b_S=',b_S)
    
    # (12/01/2020)
    # restrict to > 1
    eta = optim(eta_old, llfunc, llgrr,
                method = "L-BFGS-B", lower = 1e-6, upper = Inf)$par
    
    # fix it back... no longer fix eta to be > 0
    # eta = optim(eta_old, llfunc, llgrr, 
    #             method = "L-BFGS-B", lower = 1, upper = Inf)$par
    
    ## hack it to truth if...
    if(hack_eta){
      eta = true_Eta
    }
    
    # 2.4 if difference < tol, stop
    if(abs(beta_old - beta) < tol & abs(eta_old - eta) < tol){
      break
    }
    
  }
  
  # 3. get parameters for link activation: alphas, b_alpha
  alpha = numeric(6)
  names(alpha) = c('HH0','HI0','II0','HH1','HI1','II1')
  
  for(type in names(alpha)){
    total_dur = sum(net_c_table[[type]])
    if(!is.na(total_dur) && length(total_dur) > 0 && total_dur > 0){
      alpha[type] = C_all[type]/total_dur
    }else{
      alpha[type] = 0
    }
  }
  
  for(it in 1:maxIter){
    # 3.1 estimate b_alpha
    all_weights = c(as.matrix(net_c_table[2:7]) %*% alpha)
    has_weight = which(all_weights > 0)
    XX_alpha = XX[has_weight,]
    all_weights = all_weights[has_weight]
    alpha_poi = glm(net_c_table$Nc[has_weight]~XX_alpha-1+offset(log(all_weights)), 
                    family = poisson())
    b_alpha = alpha_poi$coefficients %>% as.numeric()
    
    # hack it to zero if...
    if(hack0){
      b_alpha = rep(0, length(b_alpha))
    }
    
    # 3.2 estimate alpha
    alpha_old = alpha
    pair_effects = exp(c(XX_alpha * b_alpha))
    for(type in names(alpha)){
      total_dur = sum(net_c_table[[type]][has_weight] * pair_effects)
      if(!is.na(total_dur) && length(total_dur) > 0 && total_dur > 0){
        alpha[type] = C_all[type]/total_dur
      }else{
        alpha[type] = 0
      }
    }
    
    # 3.3 if difference < tol, stop
    if(all(abs(alpha_old - alpha) < tol)){
      break
    }
  }
  
  
  # 4. get parameters for link termination: omegas, b_omega
  omega = numeric(6)
  names(omega) = c('HH0','HI0','II0','HH1','HI1','II1')
  
  for(type in names(omega)){
    total_dur = sum(net_d_table[[type]])
    if(!is.na(total_dur) && length(total_dur) > 0 && total_dur > 0){
      omega[type] = D_all[type]/total_dur
    }else{
      omega[type] = 0
    }
  }
  
  for(it in 1:maxIter){
    # 3.1 estimate b_omega
    all_weights = c(as.matrix(net_d_table[2:7]) %*% omega)
    has_weight = which(all_weights > 0)
    XX_omega = XX[has_weight,]
    all_weights = all_weights[has_weight]
    omega_poi = glm(net_d_table$Nd[has_weight]~XX_omega-1+offset(log(all_weights)), 
                    family = poisson())
    b_omega = omega_poi$coefficients %>% as.numeric()
    
    # hack it to zero if...
    if(hack0){
      b_omega = rep(0, length(b_omega))
    }
    
    # 3.2 estimate omega
    omega_old = omega
    pair_effects = exp(c(XX_omega * b_omega))
    for(type in names(omega)){
      total_dur = sum(net_d_table[[type]][has_weight] * pair_effects)
      if(!is.na(total_dur) && length(total_dur) > 0 && total_dur > 0){
        omega[type] = D_all[type]/total_dur
      }else{
        omega[type] = 0
      }
    }
    
    # 3.3 if difference < tol, stop
    if(all(abs(omega_old - omega) < tol)){
      break
    }
  }
  
  
  return(list(beta=beta, eta=eta, b_S=b_S, gamma=gamma, 
              phi=phi, p_s = p_s,
              alpha=alpha, b_alpha=b_alpha, omega=omega, b_omega=b_omega))
}



#####
# II. for missing data

# (1) missing recovery time

# i) given exposure times etc., get the local neighborhood of person i
#    at time t^E:
#    a vector of Ia people, and a vector of Is people

# this function: need to make sure i has manifested!! (had event 9 or 10)
get_nei_expo_i <- function(i, expo_time, G_i, events, report, times){
  # expo_time: the exposure time of i
  # all adjmats at beginning points of the intervals 
  # report: the matrix of everyone's status at each report point
  # times: the report time points
  
  # the left endpoint of interval that contains expo_time
  lb = max(times[times <= expo_time])
  ix = max(which(times <= expo_time))
  
  # taken out verbose info
  #cat('lower bound for expo time', expo_time, 'is', lb,'\n')
  
  # init epid vector to record disease status
  # only need to start from the LB of the interval for new manifestation events
  #G_i = G_all[[as.character(lb)]][i,]
  N = length(G_i)
  epid = report[ix,]
  #epid[I0] = -1 # first I0: exposed, not yet infectious!
  
  # neighborhood (ith row of G, the adjmat)
  nei = G_i
  # make sure (i,i) entry is 0
  nei[i] = 0
  
  # storage vectors for Ia and Is people
  # Ia = NULL
  # Is = NULL
  
  # get events up to expo_time:
  # get all events that are related to i (both epi and net)
  # AND all the manisfestation events (event %in% c(9,10))
  # ONLY deal with events after lb (start point of the interval)
  events = events %>% 
    filter(time <= expo_time & time >= lb) %>%
    filter(per1 == i | per2 == i | event %in% c(9,10))
  
  # if there is nothing, then NO new 
  # activity related to i from lb to expo_time
  if(nrow(events) > 0){
    # go through the events one by one
    t_pre = 0
    for(r in 1:nrow(events)){
      z = events$event[r]
      t_cur = events$time[r]
      
      # then process changes to the system
      if (z %in% c(9,10)){
        # manifestation: becoming Ia or Is
        # only deal with new infections after lb
        p1 = events$per1[r]
        epid[p1] = ifelse(z==9, 1, 2)
      }else{
        # some edge stuff
        p1 = events$per1[r]
        p2 = events$per2[r]
        
        if(z %in% c(3:5)){
          # reconnection
          if(p1 == i){
            nei[p2] = 1
          }else if(p2 == i){
            nei[p1] = 1
          }
        }else if(z %in% c(6:8)){
          # disconnection
          if(p1 == i){
            nei[p2] = 0
          }else if(p2 == i){
            nei[p1] = 0
          }
        }
      }
      # set t_pre to t_cur
      t_pre = t_cur
    }
    
  }
  
  Ia = which(nei * (epid==1) == 1)
  Is = which(nei * (epid==2) == 1)
  
  list(Ia = Ia, Is = Is)
}

# ii) get the neiborhood for all people who were ever infectious
get_nei_expo_all <- function(G_all, events, expo_times, report, times){
  # expo_times: a list of
  ## exposed: ids of all people who need exposure time imputation
  ## times: currently imputed exposure times
  exposed = expo_times$exposed
  M = ifelse(length(exposed) > 100, length(exposed), 100)
  res = foreach(i=exposed, .combine = 'list', .maxcombine = M) %dopar% {
    
    expo_time = expo_times$times[exposed == i]
    
    # the left endpoint of interval that contains expo_time
    lb = max(times[times <= expo_time])
    # get i's local adjmat at lb
    G_i = G_all[[as.character(lb)]][i,]
    
    get_nei_expo_i(i, expo_time, G_i, events, report, times)
  }
  names(res) = sapply(exposed, as.character)
  res
}


# # 10/17/2020: try it out
# setwd('~/Documents/Research_and_References/Hetero_EpiNet_2020/')
# events = read_csv('hetero_ex2_dat.csv')
# names(events) = c('time','event','per1','per2')
# I0 = 103 # directly pulled from Python console
# change_times = c(5,30)
# ## also need to re-label people with 1-indexing (rather than 0-indexing)
# events$per1 = events$per1 + 1
# events$per2 = events$per2 + 1
# 
# G0 = as.matrix(read_delim('hetero_ex2_G0.txt', delim=" ", col_names = FALSE))
# colnames(G0) = NULL
# X = as.matrix(read_delim('hetero_ex2_X.txt', delim=" ", col_names = FALSE))
# colnames(X) = NULL

# get an `exposure_times` list
# exposure_times = events %>% filter(event == 1) %>% 
#   mutate(exposed = per1, times = time) %>%
#   select(exposed, times) %>% as.list()

# # try the `get_nei_expo_all` function
# infec_nei = get_nei_expo_all(G_report, events, exposure_times, 
#                              miss_dat$report, miss_dat$report.times)
# 
# # check out a few to see if their real infector is in the set
# events %>% filter(event==1, per1 %in% c(124,28,172))
# 
# # compare it with epi_tables to validate
# # sanity check: potential infector set should be larger than real infector set!
# exposed = exposure_times$exposed
# for(e in exposed){
#   Ia_ti = epi_tables[e,'Ia_ti']
#   Is_ti = epi_tables[e,'Is_ti']
#   size_a = length(infec_nei[as.character(e)]$Ia)
#   size_s = length(infec_nei[as.character(e)]$Is)
#   cat('Person',e,'Ia size:', Ia_ti>=size_a, 'Is size:', Is_ti >=size_s,'\n')
# }

# benchmark it...
# st = Sys.time()
# for(j in 1:200){
#   infec_nei = get_nei_expo_all(G0, I0, events, exposure_times)
#   cat(j,'\r')
# }
# diff = Sys.time() - st
# cat('Time taken:',diff,'\n')
# 100 times take about 1 minute... (not super bad)


# iii) function to obtain intervals to impute recovery times on, and ids of involved individuals
get_miss_recov <- function(report, times, events){
  # data storage
  lbs = NULL; ubs = NULL
  miss_recov = list()
  # go through all intervals and check
  nt = length(times)
  for(ix in 2:nt){
    lb = times[ix-1]; ub = times[ix]
    #epid.change = report[ix,] - report[ix-1,]
    #recovered = which(epid.change < 0)
    
    # 10/19/20 update: R=-2
    epid.change = report[ix,] - report[ix-1,]
    recovered = which(epid.change < 0 & report[ix,] == -2)
    
    cat("Interval:",lb,ub,"\n")
    cat("Recovered:",recovered,"\n")
    
    if(length(recovered) > 0){
      # check if they are included in the event log already
      st = min(which(events$time > lb)); en = max(which(events$time <= ub))
      events.sel = events[st:en,]
      events.sel = events[events$event == 2,]
      if(nrow(events.sel) > 0){
        exact.recovered = events.sel$per1
        recovered = recovered[! recovered %in% exact.recovered]
      }
      # record if 'recovered' non-empty
      if(length(recovered) > 0){
        lbs = c(lbs,lb); ubs = c(ubs, ub)
        miss_recov = append(miss_recov, list(recovered))
      }
    }
  }
  return(list(intervals = data.frame(lb=lbs, ub=ubs), recover = miss_recov))
}

# iv) DARCI+: lower bound update is not uniform, adapted from before
propose_recov_filter <- function(lb, ub, recovers, exposure_times, nei_infec, 
                                 gam=0.2, exp_eta=1.5){
  
  # [lb, ub]: the time interval to impute recovery times on
  # recovers: ids of people to impute recovery times for
  # nei_infec: the list of Ia and Is neighbors at time of exposure for each person
  # exposure_times: a list of
  ## exposed: ids of all people who need exposure time imputation
  ## times: currently imputed exposure times
  # gam: currect value of gamma
  # exp_eta: current value of exp(eta)
  
  # pull up infection log in this interval 
  events.infec = as.data.frame(exposure_times) %>% 
    select(per1 = exposed, time = times) %>%
    filter(time > lb & time <= ub) %>%
    arrange(time)
  
  # if events.infec non-empty:
  # obtain adjusted feasible sampling lower bounds for each person in `recovers`
  # if about-to-recover people are the only sick neighbors of someone infected at t, 
  # then randomly select one of them to mandately recover after t
  # according to the prob. in manuscript
  bounds = rep(lb, length(recovers))
  
  if(nrow(events.infec) > 0){
    for(r in 1:nrow(events.infec)){
      p1 = events.infec$per1[r]
      if(p1 %in% recovers){
        t = events.infec$time[r]
        bounds[recovers==p1] = max(bounds[recovers==p1],t)
      }
      nei_s = nei_infec[[as.character(p1)]]$Is; size_s = length(nei_s)
      nei_a = nei_infec[[as.character(p1)]]$Ia; size_a = length(nei_a)
      nei = c(nei_s, nei_a)
      poi = recovers[recovers %in% nei]
      if(length(poi)==length(nei)){
        #cat("POIs for individual",p1,":\n",poi,"\n")
        t = events.infec$time[r]
        if(size_s == 0 | size_a == 0){
          # one of Ia/Is is empty set
          poi = nei
        }else{
          # both are non-empty
          # then choose with `exp_eta` factored in
          wset = sample(1:2, 1, prob = c(size_a, size_s * exp_eta))
          if(wset == 1){
            poi = nei_a
          }else{
            poi = nei_s
          }
        }
        ## 10/19/2020: small change to maybe resolve poi==NULL case
        if(length(poi) > 1){
          if(length(poi)==1){
            p = poi
          }else{
            p = sample(poi, 1)
          }
          bounds[recovers==p] = max(bounds[recovers==p],t)
        }
      }
    }
  }
  
  # sample recovery times under the adjusted bounds
  cands = sam_texp(length(recovers), gam, bounds-lb, ub-lb) + lb
  
  # cat("Interval: [",lb,",",ub,"]\n")
  # cat("To recover:",recovers,"\n")
  # cat("Feasible lower bounds:\n",bounds,"\n")
  # cat("Proposed recovery times:\n",cands,"\n")
  
  cands
}

## try this out
# seems to work!
# LB = 28; UB = 35
# recovs = events %>% filter(time > LB & time <= UB, event==2) %>%
#   select(per1) %>% pull()
# cands = propose_recov_filter(LB, UB, recovs, exposure_times, infec_nei)


# (2) missing exposure time

# get a list of recovers - times
# recovery_times = events %>% filter(event==2) %>% 
#   select(recov=per1, time)

# i) a function to combine two data frames and sort by time
# df_combine <- function(event1, event2){
#   n1 = nrow(event1)
#   n2 = nrow(event2)
#   if(n1 == 0){
#     event2
#   }else if(n2 == 0){
#     event1
#   }else{
#     n = n1 + n2
#     df = NULL
#     i = 1; j = 1
#     times1 = event1$time
#     times2 = event2$time
#     for(k in 1:n){
#       if(times1[i] < times2[j]){
#         df = rbind(df, )
#       }
#     }
#   }
# }

# ii) for each person i, get risk function change points and values
#    on interval [expo_time-t_max, expo_time-t_min]
#    given current recovery times and exp_eta

# 11/08/2020
# changes! previously it was most likely wrong

# 04/21/2021
# hacky resolution of getting NA in the exposure time imputation...

# 08/16/2021
# track rejection sampler acceptance ratio
get_expo_risk_i <- function(i, t_i, G_all, tmax, tmin, report, times,
                            events, recovery_times, x_i, b_S = c(0,0),
                            exp_eta=1.5, beta = 0.2, phi = 0.2,
                            details=TRUE, track_accept = TRUE){
  # t_i: manifestation time of i
  # G_all: all adjmats at report times
  # recovery_times: a data frame of imputed recovery times (no need to be sorted)
  
  # x_i: i's covar vector
  # b_S: current value of b_S
  # exp_eta: current value of exp(eta)
  # beta: current value of beta
  # phi: current value of phi
  
  # details: if TRUE, return all stuff; o.w., return sampled expo time only
  
  st = max(0,t_i - tmax)
  en = ifelse(t_i - tmin <= 0, t_i, t_i - tmin)
  # the end of interval is set to t_i instead if t_i - tmin <= 0
  
  # the left endpoint of interval that contains st
  lb = max(times[times <= st])
  ix = max(which(times <= st))
  
  # init epid vector to record disease status
  # only need to start from the LB of the interval for new manifestation events
  G_i = G_all[[as.character(lb)]][i,]
  nei = G_i
  nei[i] = 0
  epid = report[ix,]
  
  #cat('length of nei:', length(nei), "length of epid", length(epid),'\n')
  
  # filter events related to i
  # AND all the manifestation events!!
  events = events %>% filter(time >= lb & time <= en) %>% 
    filter(per1 == i | per2 == i | event %in% c(9,10))
  
  # get people recovered before st
  recovered_before_st = recovery_times %>%
    filter(time >= lb & time <= st) %>%
    select(recov) %>% pull()
  
  # recovery_times = recovery_times %>% 
  #   filter(time > st & time <= en) %>%
  #   select(time=times, per1=recov)
  
  # get recovery events in (st,en]
  recovery_times = recovery_times %>% 
    filter(time > st & time <= en) %>%
    mutate(event=2, per2 = NA) %>%
    select(time, event, per1=recov, per2)
  
  # combine events with recovery_times and rank by time
  events = rbind(events, recovery_times) %>% 
    arrange(time)
  
  # get the G_i and epid at the time of st
  t = lb
  r = 1
  if(st > lb){
    while((t < st) & (r <= nrow(events))){
      t = events$time[r]
      z = events$event[r]
      p1 = events$per1[r]
      p2 = events$per2[r]
      if(t >= st){
        break
      }
      if(z %in% c(3:5)){
        # reconnection
        if(p1 == i){
          nei[p2] = 1
        }else if(p2 == i){
          nei[p1] = 1
        }
      }else if (z %in% c(6:8)){
        # disconnection
        if(p1 == i){
          nei[p2] = 0
        }else if(p2 == i){
          nei[p1] = 0
        }
      }else if (z %in% c(9,10)){
        epid[p1] = ifelse(z==9, 1, 2)
      }
      r = r+1
    }
    
    # finally account for those who recovered in [lb,st]
    epid[recovered_before_st] = -2
    
  }
  
  # now go from st to en and register all changes
  change_times = st
  risk_val = NULL
  
  # if already exhausted events
  # OR if there is no events ...
  if(r>=nrow(events)){
    change_times = c(change_times, en)
    risk_val = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta
  }else{
    prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta
    for(k in r:nrow(events)){
      t = events$time[k]
      z = events$event[k]
      p1 = events$per1[k]
      p2 = events$per2[k]
      
      if(z %in% c(3:5)){
        # reconnection
        poi = ifelse(p1==i, p2, p1)
        nei[poi] = 1
        # if this person is Is or Ia: register change point
        if(epid[poi] %in% c(1,2)){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z %in% c(6:8)){
        # disconnection
        poi = ifelse(p1==i, p2, p1)
        nei[poi] = 0
        # if this person is Is or Ia: register change point
        if(epid[poi] %in% c(1,2)){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z %in% c(9,10)){
        # manifestation
        epid[p1] = ifelse(z==9, 1, 2)
        # if p1 is connected to i
        if(nei[p1] == 1){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }else if (z == 2){
        # recovery
        epid[p1] = -2
        # if p1 is connected to i
        if(nei[p1] == 1){
          risk_val = c(risk_val, prev_risk)
          change_times = c(change_times, t)
          prev_risk = sum(nei * (epid==1)) + sum(nei * (epid==2)) * exp_eta 
        }
      }
    }
    
    risk_val = c(risk_val, prev_risk)
    change_times = c(change_times, en)
  }
  
  # get the lengths of constant-risk intervals
  # compute relative risk
  # and sample a time according to the probs
  lens = diff(change_times)
  nL = length(lens)
  if(nL == 1){
    # prev. version
    # res = runif(1, min=st, max=en)
    
    # 11/08/2020 fix
    actual_risk = risk_val * beta * exp(sum(x_i * b_S))
    res = sam_texp(1, actual_risk, a=st, b=en)
    # 08/16/2021 track acceptance/rejection of proposals
    props = 1
    # 03/10/2021: hacky resolution of NA issues
    if(is.na(res)|!is.finite(res)){res = runif(1,st,en)}
    
    # rejection sampling to scale w.r.t. latency duration
    ratio = exp(-phi * (t_i - res))
    flag = runif(1)
    while(flag > ratio){
      res = sam_texp(1, actual_risk, a=st, b=en)
      # 08/16/2021 track acceptance/rejection of proposals
      props = props + 1
      # hacky resolution of NA issues
      if(is.na(res)|!is.finite(res)){res = runif(1,st,en)}
      ratio = exp(-phi * (t_i - res))
      flag = runif(1)
    }
    
  }else{
    # previous version
    # weights = lens * risk_val
    # A = sample(nL, 1, prob = weights)
    # res = runif(1, min=change_times[A], max=change_times[A+1])
    
    # 11/08/2020 fix
    ## get the integral on each interval
    weights = lens * risk_val
    ## multiply by beta and exp(delta_i)
    ## pad a zero at beginning
    cums = c(0,cumsum(weights)) * beta * exp(sum(x_i * b_S))
    ## take the exponent of negative
    exp_cums = exp(-cums)
    ## take diff to get prob on each interval
    probs_interval = exp_cums[1:nL] - exp_cums[2:(nL+1)]
    ## sample interval, get end points of it
    A = sample(nL, 1, prob = probs_interval)
    #A_left = change_times[A]
    #A_right = change_times[A+1]
    ## sample time on that interval
    res = sam_texp(1, risk_val[A] * beta * exp(sum(x_i * b_S)), 
                   a=change_times[A], b=change_times[A+1])
    # 08/16/2021 track acceptance/rejection of proposals
    props = 1
    # 03/10/2021: hacky stuff to resolve issues of NA
    if(is.na(res)|!is.finite(res)){res = runif(1,st,en)}   
    
    # rejection sampling to scale w.r.t. latency duration
    ratio = exp(-phi * (t_i - res))
    flag = runif(1)
    while(flag > ratio){
      A = sample(nL, 1, prob = probs_interval)
      res = sam_texp(1, risk_val[A] * beta * exp(sum(x_i * b_S)), 
                     a=change_times[A], b=change_times[A+1])
      # 08/16/2021 track acceptance/rejection of proposals
      props = props + 1
      # same hacky resolution here...
      if(is.na(res)|!is.finite(res)){res = runif(1,st,en)}
      ratio = exp(-phi * (t_i - res))
      flag = runif(1)
    }
  }
  
  # return the sampled exposure time
  if(details){
    list(changes = change_times, lengths = lens, risks = risk_val, samp = res)
  }else{
    if(track_accept){
      list(res = res, props = props)
    }else{
    res
    }
  }
  
}

# # try it out
# i = 23
# t_i = events %>% filter(event %in% c(9,10), per1==i) %>%
#   select(time) %>% pull()
# er_i = get_expo_risk_i(i,t_i, G_report, tmax=7, tmin=1,
#                        miss_dat$report, miss_dat$report.times, miss_dat$events,
#                        recovery_times, details = FALSE)
# 
# i = 77 # problem with this one!! --> probably fixed!
# t_i = events %>% filter(event %in% c(9,10), per1==i) %>% 
#   select(time) %>% pull()
# er_i = get_expo_risk_i(i,t_i, G_report, tmax=7, tmin=0,
#                        miss_dat$report, miss_dat$report.times, miss_dat$events,
#                        recovery_times)
# events %>% filter(event ==1 , per1==i)
# events %>% filter(event %in% c(9,10), per1==65) # 17.8

# i = 73
# t_i = events %>% filter(event %in% c(9,10), per1==i) %>% 
#   select(time) %>% pull()
# er_i = get_expo_risk_i(i,t_i, G_report, tmax=7, tmin=0,
#                        miss_dat$report, miss_dat$report.times, miss_dat$events,
#                        recovery_times, details = TRUE)
# events.orig %>% filter(per1==i|per2==i|event %in% c(9,10))


# iii) get sampled exposure times for all people (who ever manifested)
# 08/16/2021: add tracking of acceptance ratior for rejection sampler
get_expo_times <- function(manifest_times, G_all, tmax, tmin, report, times,
                           events, recovery_times, X, b_S = c(0,0),
                           exp_eta=1.5, beta = 0.2, phi = 0.2, track_accept=TRUE){
  
  # manifest_times: list of manifest times
  ## manifested: ids of people
  ## times: time points
  # X: the covar matrix (n by p)
  
  manifested = manifest_times$manifested
  #M = ifelse(length(manifested) > 100, length(manifested), 100)
  #res = foreach(i=manifested, .combine = 'list', .maxcombine = M) %dopar% {
  # res = foreach(i=manifested, .combine = 'c') %dopar% {
  #   
  #   t_i = manifest_times$times[manifested == i]
  #   x_i = X[i,]
  #   get_expo_risk_i(i, t_i, G_all, tmax, tmin, report, times,
  #                   events, recovery_times, x_i, b_S,
  #                   exp_eta, beta, details = FALSE)
  # }
  
  # 03/10/2021: undo the parallel, do it one-by-one
  # to show the details
  res = NULL
  if(track_accept){accepts = NULL}
  for(i in manifested){
    cat('imputing expo time for',i,'...')
    t_i = manifest_times$times[manifested==i]
    x_i = X[i,]
    # 08/16/2021: include tracking option
    res_i = get_expo_risk_i(i,t_i,G_all,tmax,tmin,report,times,
                            events, recovery_times,x_i,b_S,
                            exp_eta,beta,details=FALSE,track_accept=track_accept)
    if(typeof(res_i) == 'list'){
      cat('imputed:',res_i$res,'\n')
      res = c(res,res_i$res)
      accepts = rbind(accepts, c(res_i$props, 1)) # add a vector: c(num_props, 1) for one iter record
    }else{
      cat('imputed:',res_i,'\n')
      res = c(res,res_i)
    }
    
  }
  
  # set name as person id
  names(res) = manifested
  if(track_accept){
    rownames(accepts) = manifested
    list(res = res, accepts = accepts)
  }else{
    res
  }
}

# 12/13/2020: write a simpler version of the above function 
# for easier debugging...
get_expo_times_debug <- function(manifest_times, G_all, tmax, tmin, report, times,
                                 events, recovery_times, X, b_S = c(0,0),
                                 exp_eta=1.5, beta = 0.2, phi = 0.2, expo_truth=NULL){
  
  # manifest_times: list of manifest times
  ## manifested: ids of people
  ## times: time points
  # X: the covar matrix (n by p)
  # truth: a data frame of the true exposure times (ordered by ids in manifested)
  
  
  manifested = manifest_times$manifested
  #M = ifelse(length(manifested) > 100, length(manifested), 100)
  #res = foreach(i=manifested, .combine = 'list', .maxcombine = M) %dopar% {
  
  # res = foreach(i=manifested, .combine = 'c') %dopar% {
  #   
  #   t_i = manifest_times$times[manifested == i]
  #   x_i = X[i,]
  #   get_expo_risk_i(i, t_i, G_all, tmax, tmin, report, times,
  #                   events, recovery_times, x_i, b_S,
  #                   exp_eta, beta, details = FALSE)
  # }
  
  # get expo times one by one
  # and also print truth if truth available
  res = NULL
  for(i in manifested){
    cat('Person',i,': ')
    t_i = manifest_times$times[manifested == i]
    x_i = X[i,]
    res_i = get_expo_risk_i(i, t_i, G_all, tmax, tmin, report, times,
                            events, recovery_times, x_i, b_S,
                            exp_eta, beta, details = FALSE)
    res = c(res, res_i)
    if(!is.null(expo_truth)){
      truth_i = expo_truth$time[expo_truth$per1 == i]
      cat('imputed exposure time=',res_i, '; truth =', truth_i,'\n')
    }else{
      cat('imputed exposure time=',res_i, '\n')
    }
    
  }
  # set name as person id
  names(res) = manifested
  res
}


## try it out
# ### get a manifested list first
# manifest = events %>% filter(event %in% c(9,10)) %>%
#   select(manifested = per1, times = time) %>%
#   as.list()
# 
# ### apply the function
# imp_expo_times = get_expo_times(manifest, G_report, 7, 0, 
#                                 miss_dat$report, miss_dat$report.times, 
#                                 miss_dat$events, recovery_times, exp_eta=1.5) 
# imp_expo_times


# (3) combine observed data with imputed recovery and exposure times
combine_data <- function(events, expo_times, recov_times){
  # events: the observed data; without exposure (1) and recovery (2) events
  # expo_times: a list of imputed exposure times
  ## exposed: id
  ## times: time
  # recov_times: a data frame of imputed recovery times (currently):
  ##  recov: id
  ##  time: time
  
  expo_times = data.frame(expo_times) %>% 
    mutate(event = 1, per2=NA) %>%
    select(time = times, event, per1=exposed, per2)
  
  recov_times = recov_times %>%
    mutate(event = 2, per2=NA) %>%
    select(time = time, event, per1=recov, per2)
  
  events_aug = rbind(events, expo_times, recov_times) %>%
    arrange(time)
  
  events_aug
}

## try it out
# events_aug = combine_data(miss_dat$events, exposure_times, recovery_times_old)



# 03/30/2021
# (1) helper functions to calculate 1st and 2nd derivatives for exp_eta
# (2) upgraded function for est_var

## exp_eta derivative; the first term only!
get_eta_1 <- function(Ia_ti, Is_ti, exp_eta){
  # Ia_ti, Is_ti: length-n vector of connected infectives at expo time
  # estimate of exp(eta) to plug in
  
  has_expo = which(Ia_ti + Is_ti > 0)
  if(length(has_expo)>0){
    sum(Is_ti[has_expo]/(Ia_ti[has_expo] + Is_ti[has_expo] * exp_eta))
  }else{
    0
  }
}

## exp_eta 2nd-order derivative (get the negative of it)
get_eta_2 <- function(Ia_ti, Is_ti, exp_eta){
  # Ia_ti, Is_ti: length-n vector of connected infectives at expo time
  # estimate of exp(eta) to plug in
  
  has_expo = which(Ia_ti + Is_ti > 0)
  if(length(has_expo)>0){
    sum(Is_ti[has_expo]^2/(Ia_ti[has_expo] + Is_ti[has_expo] * exp_eta)^2)
  }else{
    0
  }
}


# the upgraded version of "est_var"

## 04/21/2021: the simulation data version (no xi and b_E)
est_var <- function(res, est, last=20){
  ## res: list of raw results from a stochastic EM run
  ## est: list of parameter estimates (e.g., obtained from averaging)
  
  X = res$X
  Expo_Ias = res$Expo_Ias
  Expo_Iss = res$Expo_Iss
  Ia_tis = res$Ia_tis
  Is_tis = res$Is_tis
  sum_latent = res$sum_latent
  sum_sick = res$sum_sick
  counts = res$counts
  
  ## only take the last "last" samples for the summary statistics
  nsamps = nrow(Expo_Ias)
  Expo_Ias = Expo_Ias[max(nsamps-last+1,0):nsamps,]
  Expo_Iss = Expo_Iss[max(nsamps-last+1,0):nsamps,]
  Ia_tis = Ia_tis[max(nsamps-last+1,0):nsamps,]
  Is_tis = Is_tis[max(nsamps-last+1,0):nsamps,]
  sum_latent = sum_latent[max(nsamps-last+1,0):nsamps]
  sum_sick = sum_sick[max(nsamps-last+1,0):nsamps]
  
  ## (1) beta, exp_eta, and b_S
  beta = est$beta
  exp_eta = est$exp_eta
  b_S = est$b_S
  
  indiv_effects = exp(c(X %*% b_S))
  
  ### the b_S part
  grads_b_S = -beta * (Expo_Ias + exp_eta * Expo_Iss) %*% diag(indiv_effects) %*% X
  
  mean_expo_Ia = apply(Expo_Ias, 2, mean, na.rm=TRUE)
  mean_expo_Is = apply(Expo_Iss, 2, mean, na.rm=TRUE)
  indiv_expos = indiv_effects * (mean_expo_Ia + mean_expo_Is * exp_eta)
  
  negHess_b_S = beta * t(X) %*% diag(indiv_expos) %*% X
  
  ### 1st-order der. of beta and exp_eta
  sum_Is = (Expo_Iss %*% indiv_effects) %>% as.vector()
  sum_Ia = (Expo_Ias %*% indiv_effects) %>% as.vector()
  
  grads_beta = counts['nE']/beta - (sum_Ia + exp_eta * sum_Is)
  grads_exp_eta = as.vector(sapply(1:nrow(Expo_Iss), 
                                   function(j) get_eta_1(Ia_tis[j,], Is_tis[j,], exp_eta))) -
    beta * sum_Is
  
  gradients = cbind(grads_beta, grads_exp_eta, grads_b_S)
  
  ### the cross-terms for Hessian (also get the negative of it...)
  cross_beta_eta = sum_Is
  cross_eta_b_S = (beta * t(X) %*% (indiv_effects * mean_expo_Is)) %>% as.vector()
  cross_beta_b_S = (t(X) %*% indiv_expos) %>% as.vector()
  
  ### the (negative of) 2nd order terms for beta and exp_eta
  negHess_beta = counts['nE']/beta^2
  negHess_exp_eta = mean(sapply(1:min(nsamps,last), 
                                function(j) get_eta_2(Ia_tis[j,], Is_tis[j,], exp_eta)),
                         na.rm = TRUE)
  
  ### build the (negative of) mean Hessian matrix for (beta, exp_eta, b_S)
  p = ncol(X)
  negHess = matrix(0, p+2, p+2)
  negHess[3:(p+2),3:(p+2)] = negHess_b_S
  negHess[1:2,1:2] = matrix(c(negHess_beta,cross_beta_eta,
                              cross_beta_eta,negHess_exp_eta),
                            2,2)
  negHess[1:2,3:(p+2)] = rbind(cross_beta_b_S,cross_eta_b_S)
  negHess[3:(p+2),1:2] = cbind(cross_beta_b_S,cross_eta_b_S)
  
  #b_S_I = negHess - cov(gradients)
  #b_S_var = diag(solve(b_S_I))
  
  expo_I = negHess - cov(gradients)
  expo_var = diag(solve(expo_I))
  
  ## (2) p_s
  p_s = est$p_s
  
  p_s_var = 1/(counts['nIs']/p_s^2 + counts['nIa']/(1-p_s)^2)
  
  ## (3) phi and gamma
  gamma = est$gamma
  phi = est$phi
  
  grads_gamma = counts['nR']/gamma - sum_sick
  #gamma_I = counts['nR']/gamma - var(grads_gamma)
  gamma_var = 1/(counts['nR']/gamma^2 - var(grads_gamma))
  
  # grads_phi = counts['nI_int']/phi - sum_latent
  # phi_var = 1/(counts['nI_int']/phi^2 - var(grads_phi))

  nI = counts['nIa'] + counts['nIs']
  grads_phi = nI/phi - sum_latent
  phi_var = 1/(nI/phi^2 - var(grads_phi))
  
  list(expo_I = expo_I, expo_var = expo_var,
       p_s_var = p_s_var,
       gamma_var = gamma_var, phi_var = phi_var)
}


# 03/30/2021: a helper function to pool results together across multiple runs

# 04/29/2021: slight modification to include averaged estimates as well
# (only average over last iters across the chains)
# (and ONLY focus on the 1-d params)

# 08/05/2021: 
# slight change to record the diff from truth for b_S
pool_res <- function(res_dir, r_ls, last=20, diff_bS = TRUE){
  Expo_Ias = NULL
  Expo_Iss = NULL
  Ia_tis = NULL
  Is_tis = NULL
  sum_latent = NULL
  sum_sick = NULL
  
  avg_est = list()
  vnames = c("beta", "eta", "exp_eta", "gamma", "phi", "p_s")
  for(v in vnames){avg_est[[v]] = 0}
  
  ## add b_S
  avg_est$b_S = rep(0,2)
  
  for(r in r_ls){
    #fname = paste0(res_dir, r,'.rds')
    fname = paste0(res_dir, r)
    res = readRDS(fname)
    nsamps = nrow(res$Expo_Ias)
    Expo_Ias = rbind(Expo_Ias, res$Expo_Ias[max(nsamps-last+1,0):nsamps,])
    Expo_Iss = rbind(Expo_Iss, res$Expo_Iss[max(nsamps-last+1,0):nsamps,])
    Ia_tis = rbind(Ia_tis, res$Ia_tis[max(nsamps-last+1,0):nsamps,])
    Is_tis = rbind(Is_tis, res$Is_tis[max(nsamps-last+1,0):nsamps,])
    sum_latent = c(sum_latent, res$sum_latent[max(nsamps-last+1,0):nsamps])
    sum_sick = c(sum_sick, res$sum_sick[max(nsamps-last+1,0):nsamps])
    
    ## get last iter as well
    for(v in vnames){
      avg_est[[v]] = avg_est[[v]] + res$params[[v]][nsamps]
    }
    ## also do b_S
    if(diff_bS){
      ## get b_S truth first
      b_S_tru = res$truth$b_S
      ## and calculate difference from the truth
      this_b_S = res$params$b_S[nsamps,] - b_S_tru
    }else{
      this_b_S = res$params$b_S[nsamps,]
    }
    
    avg_est$b_S = avg_est$b_S + this_b_S
  }
  
  # take average over the last iter estimates
  for(v in vnames){avg_est[[v]] = avg_est[[v]]/length(r_ls)}
  avg_est$b_S = avg_est$b_S/length(r_ls)
  
  list(Expo_Ias = Expo_Ias, Expo_Iss = Expo_Iss, 
       Ia_tis = Ia_tis, Is_tis = Is_tis, 
       sum_latent = sum_latent, sum_sick = sum_sick,
       counts = res$counts, X = res$X,
       estimates = avg_est)
}
