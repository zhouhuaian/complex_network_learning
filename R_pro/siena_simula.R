library('RSiena')
library(sna)
library(plyr)
library(Matrix)
options(stringsAsFactors=F)

siena07ToConvergence <- function(alg, dat, eff, nbrNodes, prev=NULL){
  numr <- 0
  myeff <- eff
  if(!is.null(prev)) {
    myeff <- updateTheta(myeff, prev)
  }
  # Score-type tests go here.
  # e.g.
  # myeff <- setEffect(myeff, inPopSqrt, fix=T, test=T, initialValue=0)
  
  ans <- siena07(alg, data=dat, effects=myeff, nbrNodes=nbrNodes, verbose=F, useCluster=T, returnDeps=F) # the first run
  repeat {
    numr <- numr+1
    # count number of repeated runs
    maxt <- max(abs(ans$tconv[!myeff$fix[myeff$include]]))
    # convergence indicator, excluding the fixed effects
    cat(numr, maxt,"\n")
    # report how far we are
    if (maxt < 0.10) {break} # success
    if (maxt > 5) {break}
    # divergence without much hope
    # of returning to good parameter values
    if (numr > 20) {break} # now it has lasted too long

    # Need score-type tests both above and below
    # e.g.
    # myeff <- setEffect(myeff, inPopSqrt, fix=T, test=T, initialValue=0)
    
    ans <- siena07(alg, data=dat, effects=myeff, nbrNodes=nbrNodes, prevAns=ans, verbose=T, useCluster=T, returnDeps=F)
  }
  print(numr)
  return(ans)
}

# Find collinear parameters
# from ./cov.dat covariance
# matrix and a given threshold
find_coll <- function(ans, thresh) {
  cov_mat <- as.matrix(read.table('./cov.dat'))
  cov_tril <- tril(cov_mat, k=-1)
  ind <- get_ind(cov_tril, thresh)
  apply(ind, 1, function(row) c(ans$effects$effectName[row[1]], ans$effects$effectName[row[2]]))
}

GeodesicDistribution <- function (i, data, sims, period, groupName,
                                  varName, levls=c(1:5,Inf), cumulative=TRUE, ...) {
  x <- networkExtraction(i, data, sims, period, groupName, varName)
  require(sna)
  a <- sna::geodist(symmetrize(x))$gdist
  if (cumulative)
  {
    gdi <- sapply(levls, function(i){ sum(a<=i) })
  }
  else
  {
    gdi <- sapply(levls, function(i){ sum(a==i) })
  }
  names(gdi) <- as.character(levls)
  gdi
}

# Helper function
get_ind <- function(mat, thresh) {
  ret <- c()
  for(i in 1:nrow(mat)) {
    for(j in 1:ncol(mat)) {
      if(!is.na(mat[i,j]) & abs(mat[i,j]) > thresh & mat[i,j] != 0) {
        print(mat[i,j])
        print(i)
        print(j)
        ret <- rbind(ret, c(i,j))
      }
    }
  }
  return(ret)
}

retrieveDataMatrix <- function(data_name, remove_waves, remove_devs) {
  mat <- as.matrix(read.table(construct_dat_file_path(project_name, data_name)))
  if(!is.null(remove_waves)) {
    mat <- mat[,-remove_waves]
  }
  if(!is.null(remove_devs)) {
    mat <- mat[-remove_devs,]
  }
  mat
}

# Call BEFORE remove_waves is removed from num_waves
createTimeDummyMatrix <- function(nwaves, ndevs, remove_waves_listing, type) {
  mat <- matrix(0, nrow=ndevs, ncol=nwaves)
  for(i in 1:length(remove_waves_listing)) {
    rwave <- remove_waves_listing[[i]]$wave
    if(rwave - 1 < 1) {
      next
    }
    if(remove_waves_listing[[i]]$type == type) {
      mat[,rwave-1] <- mat[,rwave-1] + 1
    }
    mat <- mat[,-rwave]
    for(j in 1:length(remove_waves_listing)) {
      if(remove_waves_listing[[j]]$wave > rwave) {
        remove_waves_listing[[j]]$wave <- remove_waves_listing[[j]]$wave - 1
      }
    }
  }
  mat
}

# all_data is MUST BE ALREADY FILTERED
# for the proper number of waves.
# There isn't a good workaround for this right now.
processBehNormalizedDataWithLeavers <- function(all_data, first_commit_data, last_commit_data, ndevs, behavior_name, remove_devs) {
  dates <- unique(all_data$ownership_date)
  beh_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  names(last_commit_data) <- c("committer_id", "max")
  date_ind <- 0
  for(dat in dates) {
    date_ind <- date_ind + 1
    print(date_ind)
    committer_ids <- all_data[all_data$ownership_date == dat,]$committer_id
    total_df <- all_data[all_data$ownership_date == dat,c("ownership_date","committer_id",behavior_name)]
    
    norm_val <- data.frame(norm_val=as.numeric(as.POSIXct(dat) - as.POSIXct(first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"min"])),
                           committer_id=first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"committer_contributor_id"])
    if(any(norm_val < 0)) {
      print("NORMALIZED VALUE < 0. ABORTING.")
      return
    }
    total_df <- merge(total_df, norm_val)
    total_df$normalized <- total_df[[behavior_name]] / total_df$norm_val
    # Make sure to initialize the behaviors to 0 so that we can
    # then get the right indices for not including normalized values of 0.
    # i.e. the quantiles should be calculated without 0s in the normalized values
    total_df$beh <- 0
    # If < 5 quantiles, introduce a tiny bit of jitter so we can get quantiles with behavior
    # values of the same scale as the actual normalized values
    if(length(unique(quantile(total_df[total_df$normalized > 0, ]$normalized, probs=seq(0,1,by=.25)))) < 5) {
      total_df$beh <- with(total_df, 
                           cut(normalized,
                               breaks=quantile(normalized + seq_along(normalized) * .Machine$double.eps, probs=seq(0,1,by=.25)), include.lowest=T))
    } else {
      total_df[total_df$normalized > 0, "beh"] <- with(total_df[total_df$normalized > 0, ], 
                                                       cut(normalized,
                                                           breaks=quantile(normalized, probs=seq(0,1,by=.25)), include.lowest=T))
    }
    total_df$beh <- as.numeric(total_df$beh)
    
    total_df <- merge(total_df, last_commit_data)
    # Save the behaviors that are going to be set to 0 in case we need them
    # in the if statement below
    saved_beh_zero <- total_df[as.POSIXct(dat) > as.POSIXct(total_df$max), "beh"]
    total_df[as.POSIXct(dat) > as.POSIXct(total_df$max), "beh"] <- 0
    
    # Now we have to recalculate the quantiles to exclude behaviors of leavers.
    # If there's < 5 quantiles, set behaviors to what they were before we
    # reset peoples' behaviors to 0.
    # We can't do the jitter approach here because we may have < 5 people that have behavior values > 0,
    # which can't be handled by the jitter approach as done above.
    # This is an edge case that has to be dealt with
    if(length(unique(quantile(total_df[total_df$beh != 0,]$normalized, probs=seq(0, 1, by=.25)))) < 5) {
      # Use this to reset the values to the values before getting rid of zeros, if that's
      # what we want to do instead
      total_df[total_df$beh == 0, "beh"] <- saved_beh_zero
    } else {
      total_df[total_df$beh != 0, "beh"] <- with(total_df[total_df$beh != 0,],
                                                 cut(normalized,
                                                     breaks=quantile(normalized, probs=seq(0,1,by=.25)), include.lowest=T))
    }
    total_df$beh <- as.numeric(total_df$beh)
    
    beh_mat[total_df$committer_id, date_ind] <- total_df$beh
  }
  if(!is.null(remove_devs)) {
    beh_mat <- beh_mat[-remove_devs,]
  }
  beh_mat
}

processBehNotNormalizedDataWithLeavers <- function(all_data, first_commit_data, last_commit_data, ndevs, behavior_name, remove_devs) {
  dates <- unique(all_data$ownership_date)
  beh_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  names(last_commit_data) <- c("committer_id", "max")
  date_ind <- 0
  for(dat in dates) {
    date_ind <- date_ind + 1
    print(date_ind)
    committer_ids <- all_data[all_data$ownership_date == dat,]$committer_id
    total_df <- all_data[all_data$ownership_date == dat,c("ownership_date","committer_id",behavior_name)]
    names(total_df) <- c("ownership_date","committer_id","behavior")
    
    # Make sure to initialize the behaviors to 0 so that we can
    # then get the right indices for not including normalized values of 0.
    # i.e. the quantiles should be calculated without 0s in the normalized values
    total_df$beh <- 0
    # If < 5 quantiles, introduce a tiny bit of jitter so we can get quantiles with behavior
    # values of the same scale as the actual normalized values
    #total_df$normalized
    if(length(unique(quantile(total_df[total_df[["behavior"]] > 0, "behavior"], probs=seq(0,1,by=.25)))) < 5) {
      total_df$beh <- with(total_df, 
                           cut(behavior,
                               breaks=quantile(behavior + seq_along(behavior) * .Machine$double.eps, probs=seq(0,1,by=.25)), include.lowest=T))
    } else {
      total_df[total_df[["behavior"]] > 0, "beh"] <- with(total_df[total_df[["behavior"]] > 0, ], 
                                                          cut(behavior,
                                                              breaks=quantile(behavior, probs=seq(0,1,by=.25)), include.lowest=T))
    }
    total_df$beh <- as.numeric(total_df$beh)
    
    total_df <- merge(total_df, last_commit_data)
    # Save the behaviors that are going to be set to 0 in case we need them
    # in the if statement below
    saved_beh_zero <- total_df[as.POSIXct(dat) > as.POSIXct(total_df$max), "beh"]
    total_df[as.POSIXct(dat) > as.POSIXct(total_df$max), "beh"] <- 0
    
    # Now we have to recalculate the quantiles to exclude behaviors of leavers.
    # If there's < 5 quantiles, set behaviors to what they were before we
    # reset peoples' behaviors to 0.
    # We can't do the jitter approach here because we may have < 5 people that have behavior values > 0,
    # which can't be handled by the jitter approach as done above.
    # This is an edge case that has to be dealt with
    if(length(unique(quantile(total_df[total_df$beh != 0,"behavior"], probs=seq(0, 1, by=.25)))) < 5) {
      # Use this to reset the values to the values before getting rid of zeros, if that's
      # what we want to do instead
      total_df[total_df$beh == 0, "beh"] <- saved_beh_zero
    } else {
      total_df[total_df$beh != 0, "beh"] <- with(total_df[total_df$beh != 0,],
                                                 cut(behavior,
                                                     breaks=quantile(behavior, probs=seq(0,1,by=.25)), include.lowest=T))
    }
    total_df$beh <- as.numeric(total_df$beh)
    
    beh_mat[total_df$committer_id, date_ind] <- total_df$beh
  }
  if(!is.null(remove_devs)) {
    beh_mat <- beh_mat[-remove_devs,]
  }
  beh_mat
}

# Used to normalize data by time
processBehNormalizedDataNoLeavers <- function(all_data, first_commit_data, last_commit_data, ndevs, behavior_name, remove_devs) {
  dates <- unique(all_data$ownership_date)
  beh_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  names(last_commit_data) <- c("committer_id", "max")
  date_ind <- 0
  for(date in dates) {
    date_ind <- date_ind + 1
    print(date_ind)
    committer_ids <- all_data[all_data$ownership_date == date,]$committer_id
    total_df <- all_data[all_data$ownership_date == date,c("ownership_date","committer_id",behavior_name)]
    
    norm_val <- data.frame(norm_val=as.numeric(as.POSIXct(date) - as.POSIXct(first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"min"])),
                           committer_id=first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"committer_contributor_id"])
    if(any(norm_val < 0)) {
      print("NORMALIZED VALUE < 0. ABORTING.")
      return
    }
    total_df <- merge(total_df, norm_val)
    total_df$normalized <- total_df[[behavior_name]] / total_df$norm_val
    total_df$beh <- with(total_df, 
                         cut(normalized,
                             breaks=quantile(normalized, probs=seq(0,1,by=.25)), include.lowest=T))
    total_df$beh <- as.numeric(total_df$beh)
    
    beh_mat[total_df$committer_id, date_ind] <- total_df$beh
  }
  if(!is.null(remove_devs)) {
    beh_mat <- beh_mat[-remove_devs,]
  }
  beh_mat
}

retrieveDataMatrixFromAllData <- function(all_data, ndevs, col_name, remove_devs) {
  dates <- unique(all_data$ownership_date)
  col_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  date_ind <- 0
  for(date in dates) {
    date_ind <- date_ind + 1
    committer_ids <- all_data[all_data$ownership_date == date,]$committer_id
    total_df <- all_data[all_data$ownership_date == date,c("ownership_date","committer_id",col_name)]
    col_mat[total_df$committer_id, date_ind] <- total_df[[col_name]]
  }
  if(!is.null(remove_devs)) {
    col_mat <- col_mat[-remove_devs,]
  }
  col_mat
}

retrieveAgeMatrix <- function(all_data, first_commit_data, ndevs, remove_devs) {
  dates <- unique(all_data$ownership_date)
  col_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  date_ind <- 0
  for(date in dates) {
    date_ind <- date_ind + 1
    committer_ids <- all_data[all_data$ownership_date == date,]$committer_id
    total_df <- all_data[all_data$ownership_date == date,c("ownership_date","committer_id")]
    
    norm_val <- data.frame(norm_val=as.numeric(as.POSIXct(date) - as.POSIXct(first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"min"])),
                           committer_id=first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"committer_contributor_id"])
    if(any(norm_val < 0)) {
      print("NORMALIZED VALUE < 0. ABORTING.")
      return
    }
    total_df <- merge(total_df, norm_val)
    
    col_mat[total_df$committer_id, date_ind] <- total_df[["norm_val"]]
  }
  if(!is.null(remove_devs)) {
    col_mat <- col_mat[-remove_devs,]
  }
  col_mat
}

retrieveAgeMatrixAllOnes <- function(all_data, first_commit_data, ndevs, remove_devs) {
  dates <- unique(all_data$ownership_date)
  col_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  date_ind <- 0
  for(date in dates) {
    date_ind <- date_ind + 1
    committer_ids <- all_data[all_data$ownership_date == date,]$committer_id
    total_df <- all_data[all_data$ownership_date == date,c("ownership_date","committer_id")]
    
    norm_val <- data.frame(norm_val=as.numeric(as.POSIXct(date) - as.POSIXct(first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"min"])),
                           committer_id=first_commit_data[first_commit_data$committer_contributor_id %in% committer_ids,"committer_contributor_id"])
    if(any(norm_val < 0)) {
      print("NORMALIZED VALUE < 0. ABORTING.")
      return
    }
    total_df <- merge(total_df, norm_val)
    
    col_mat[total_df$committer_id, date_ind] <- 1
  }
  if(!is.null(remove_devs)) {
    col_mat <- col_mat[-remove_devs,]
  }
  col_mat
}

getEmailCumulativeCountData <- function(dates, remove_devs, ndevs) {
  email_mat <- matrix(0, nrow=ndevs, ncol=length(dates))
  date_ind <- 0
  for(j in 1:length(dates)) {
    date <- dates[j]
    date_ind <- date_ind + 1
    for(i in 1:ndevs) {
      email_mat[i, date_ind] <- length(emails[as.POSIXct(emails$date) <= as.POSIXct(date) & emails$from == i,"date"])
    }
  }
  if(!is.null(remove_devs)) {
    email_mat <- email_mat[-remove_devs,]
  }
  email_mat
}

createCompositionChange <- function(dates, split_k, remove_devs, ndevs) {
  first_date <- as.POSIXct(dates[1])
  last_date <- as.POSIXct(dates[length(dates)])
  
  emails_from_max <- aggregate(as.POSIXct(emails$date), list(emails$from), max)
  
  emails_to_max <- aggregate(as.POSIXct(emails$date), list(emails$to), max)
  
  emails_max <- rbind(emails_from_max, emails_to_max)
  emails_max <- aggregate(as.POSIXct(emails_max[,2]), list(emails_max[,1]), max)
  names(emails_max) <- c("committer_id", "date")
  
  emails_from_min <- aggregate(as.POSIXct(emails$date), list(emails$from), min)
  
  emails_to_min <- aggregate(as.POSIXct(emails$date), list(emails$to), min)
  
  emails_min <- rbind(emails_from_min, emails_to_min)
  emails_min <- aggregate(as.POSIXct(emails_min[,2]), list(emails_min[,1]), min)
  names(emails_min) <- c("committer_id", "date")
  
  rowSu <- lapply(savenets, function(x) { rowSums(x) })
  colSu <- lapply(savenets, function(x) { colSums(x) })
  totSu <- vector("list", split_k)
  for(i in 1:split_k) {
    totSu[[i]] <- rowSu[[i]] + colSu[[i]]
  }
  
  totSuMat <- matrix(unlist(totSu), nrow=split_k, ncol=ndevs, byrow=T)
  
  leaves <- apply(totSuMat, 2, function(x) {
    m_leaves <- c()
    prev <- 1
    for(i in 2:length(x)) {
      if(x[prev] > 0 & x[i] == 0) {
        m_leaves <- c(m_leaves, i - .1)
      }
      prev <- i
    }
    return(m_leaves)
  })
  
  joins <- apply(totSuMat, 2, function(x) {
    m_joins <- c()
    prev <- 1
    for(i in 2:length(x)) {
      if(x[prev] == 0 & x[i] > 0) {
        m_joins <- c(m_joins, i - .9)
      }
      prev <- i
    }
    return(m_joins)
  })
  for(i in 1:ndevs) {
    if(totSuMat[1,i] > 0) {
      if(is.null(joins[[i]])) {
        joins[[i]] <- 1
      } else {
        joins[[i]] <- c(joins[[i]], 1)
      }
    }
    
    if(totSuMat[split_k,i] > 0) {
      if(is.null(leaves[[i]])) {
        leaves[[i]] <- split_k
      } else {
        leaves[[i]] <- c(leaves[[i]], split_k)
      }
    }
    if(is.null(joins[[i]])) {
      joins[[i]] <- 1
    }
    if(is.null(leaves[[i]])) {
      leaves[[i]] <- split_k
    }
  }
  comp <- rep(list(c(1, split_k)), ndevs)
  for(i in 1:ndevs) {
    comp[[i]] <- sort(c(joins[[i]], leaves[[i]]))
  }
  comp[remove_devs] <- NULL
  comp
}

getSplitDates <- function(k) {
  split_ind <- ceiling(seq(nrow(emails) / k, nrow(emails), by=nrow(emails)/k))
  split_dates <- emails[split_ind, "date"]
  return(split_dates)
}

getSplitIndsFullNet <- function(k) {
  split_dates <- getSplitDates(k)
  start_date <- min(as.POSIXct(emails$date))
  prev_date <- start_date
  end_date <- max(as.POSIXct(emails$date))
  ndays <- as.numeric(ceiling(end_date - start_date))
  one_unit <- 1*60*60*24
  full_net_dates <- c()
  for(i in 1:ndays) {
    cur_date <- prev_date + one_unit
    full_net_dates <- c(full_net_dates, cur_date)
    prev_date <- cur_date
  }
  closest <- c()
  for(i in 1:length(split_dates)) {
    chk <- as.numeric(split_dates[i])
    closest <- c(closest, which(abs(full_net_dates-chk)==min(abs(full_net_dates-chk))))
  }
  return(closest)
}

generate_commit_window_data <- function(wave_dates, rate_den, commits, ndevs, remove_devs, ret_str) {
  # window size = half life of a non-refreshed edge
  one_unit <- 1*60*60*24
  commits_in_window <- matrix(0, nrow=ndevs, ncol=length(wave_dates))
  window_size <- rate_den * log(2)
  date_ind <- 0
  for(date_ind in 1:length(wave_dates)) {
    date <- wave_dates[date_ind]
    for(dev in 1:ndevs) {
      cs <- commits[commits$committer_id == dev & as.POSIXct(commits$committer_dt) > as.POSIXct(date) - window_size*one_unit & as.POSIXct(commits$committer_dt) <= as.POSIXct(date),]
      commits_in_window[dev,date_ind] <- length(unique(cs$committer_dt))
    }
  }
  commits_in_window_df <- as.data.frame(commits_in_window)
  # I'm so sorry
  for(date_ind in 1:length(wave_dates)) {
    print(date_ind)
    whi <- which(commits_in_window_df[,date_ind] > 0)
    if(length(unique(quantile(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0], probs=seq(0,1,by=.25)))) < 5) {
      commits_in_window_df[whi,date_ind] <- cut(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0] + seq_along(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0])*.Machine$double.eps,
                                                breaks=quantile(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0] + seq_along(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0])*.Machine$double.eps, probs=seq(0,1,by=.25)),
                                                include.lowest=T)
    } else {
      commits_in_window_df[whi,date_ind] <- cut(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0],
                                                breaks=quantile(commits_in_window_df[,date_ind][commits_in_window_df[,date_ind] > 0], probs=seq(0,1,by=.25)),
                                                include.lowest=T)
    }
  }
  if(!is.null(remove_devs)) {
    commits_in_window <- commits_in_window[-remove_devs,]
    commits_in_window_df <- commits_in_window_df[-remove_devs,]
  }
  if(ret_str == "mat") {
    return(commits_in_window)
  } else if(ret_str == "beh") {
    return(commits_in_window_df)
  } else {
    stop("Invalid ret str in generate_commit_window_data")
  }
}

# Input project name here, e.g. ant
project_name <- "PROJECT_NAME_PLACEHOLDER"
# Decay rate
# MUST BE AN INTEGER
rate_den <- "DECAY_RATE_INTEGER_PLACEHOLDER"
# k is 8 for our publication - number of splits in network data, as defined in
# build_network_with_decay_function.R
k <- 8
model_name <- "MODEL_NAME_PLACEHOLDER"
# NOTE: IMPORTANT: SET THIS WHEN YOU MANUALLY KNOW WHICH WAVES WE WANT TO MODEL
# For axis2_java, this was waves 2:6
selected_waves <- 2:6

# What behavior are we looking at?
# e.g. ownership_count_70, commit_window
# Name of column from all_data
beh_name <- "ownership_count_70"

emails_file <- "EMAIL_DATA_FILE_PLACEHOLDER"
emails <- read.csv(emails_file, header=F, stringsAsFactors=F)
names(emails) <- c("from","to","date")
emails$date <- as.POSIXct(emails$date)

# getSplitDates accesses the email df globally
wave_dates <- getSplitDates(k)

commits_file <- "COMMITS_DATA_FILE_PLACEHOLDER"
commits <- read.csv(commits_file, header=T, stringsAsFactors=F)
# Convert the dates to POSIXct AFTER WE CREATE first_commit_data, last_commit_data

all_data_file <- "ALL_DATA_FILE_PLACEHOLDER"
all_data <- read.csv(file=all_data_file, sep=",", header=T)

# Load in "savenets" from
# our decay script
# build_network_with_decay_function.R
load(file="SAVENETS_DECAY_FILE_PLACEHOLDER")

num_waves <- 1:k
num_devs <- 1:nrow(savenets[[1]])
ndevs <- length(num_devs)

first_commit_data <<- data.frame(committer_contributor_id=numeric(length(unique(commits$committer_id))), min=character(length(unique(commits$committer_id))))
last_commit_data <<- data.frame(committer_contributor_id=numeric(length(unique(commits$committer_id))), max=character(length(unique(commits$committer_id))))


sapply(unique(commits$committer_id), function(dev, commits) {
  first_commit_data[dev,] <<- c(dev, min(commits[commits$committer_id == dev,"committer_dt"]))
}, commits=commits)

sapply(unique(commits$committer_id), function(dev, commits) {
  last_commit_data[dev,] <<- c(dev, max(commits[commits$committer_id == dev,"committer_dt"]))
}, commits=commits)

commits$committer_dt <- as.POSIXct(commits$committer_dt)

remove_devs <- c()
# Remove devs with no emails
remove_devs_email <- num_devs[!(num_devs %in% unique(emails$from) | num_devs %in% unique(emails$to))]
# Remove devs outside the range of email wave dates
email_range <- c(wave_dates[1], wave_dates[k])
emails_from_max <- aggregate(as.POSIXct(emails$date), list(emails$from), max)

emails_to_max <- aggregate(as.POSIXct(emails$date), list(emails$to), max)

emails_max <- rbind(emails_from_max, emails_to_max)
emails_max <- aggregate(as.POSIXct(emails_max[,2]), list(emails_max[,1]), max)
names(emails_max) <- c("committer_id", "edate")

emails_from_min <- aggregate(as.POSIXct(emails$date), list(emails$from), min)

emails_to_min <- aggregate(as.POSIXct(emails$date), list(emails$to), min)

emails_min <- rbind(emails_from_min, emails_to_min)
emails_min <- aggregate(as.POSIXct(emails_min[,2]), list(emails_min[,1]), min)
names(emails_min) <- c("committer_id", "sdate")

ranges <- join(emails_min, emails_max)
for(i in 1:ndevs) {
  m_range <- ranges[ranges$committer_id == i,c("sdate", "edate")]
  if(i %in% remove_devs_email 
     | (m_range[1,"sdate"] < email_range[1] & m_range[1, "edate"] < email_range[1]) 
     | (m_range[1,"sdate"] > email_range[2] & m_range[1, "edate"] > email_range[2])) {
    remove_devs_email <- c(remove_devs_email, i)
  }
}
# Change the behavior according to the run
# See beh_name above
if(beh_name != "commit_window") {
  beh <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, beh_name, remove_devs)
} else {
  beh <- generate_commit_window_data(wave_dates, rate_den, commits, ndevs, remove_devs, "beh")
}

# Remove devs who enter and leave the network more than once
comp <- createCompositionChange(wave_dates, k, c(), ndevs)
for(i in 1:ndevs) {
  if(length(comp[[i]]) > 2) {
    remove_devs <- c(remove_devs, i)
  }
  # also remove people who leave before the start
  # of our selected waves or enter after the
  # end of our selected waves.
  if(comp[[i]][2] < selected_waves[1] | comp[[i]][1] > selected_waves[length(selected_waves)]) {
    remove_devs <- c(remove_devs, i)
  }
}
remove_devs <- sort(unique(remove_devs))


degchk <- list()
for(i in 1:k) {
  degchk[[i]] <- rowSums(savenets[[i]]) + colSums(savenets[[i]])
}

complete_degs <- matrix(unlist(degchk), nrow=k, ncol=length(num_devs))
remove_devs <- c(remove_devs, which(colSums(complete_degs) == 0))
remove_devs <- sort(unique(remove_devs))

max_outdegs <- list()
for(i in 1:k) {
  max_outdegs[[i]] <- max(rowSums(savenets[[i]]))
}

max_indegs <- list()
for(i in 1:k) {
  max_indegs[[i]] <- max(colSums(savenets[[i]]))
}

if(length(remove_devs) == 0) {
  remove_devs <- NULL
}

if(!is.null(remove_devs)) {
  num_devs <- num_devs[-remove_devs]
}

network_vector <- c()
for(i in num_waves) {
  net <- savenets[[i]]
  net <- net[num_devs, num_devs]
  network_vector <- c(network_vector, net)
}

email_net_array <- array(network_vector, dim=c(length(num_devs), length(num_devs), length(num_waves)))
email_network <- sienaDependent(email_net_array[,,selected_waves], allowOnly=T)

comp <- createCompositionChange(wave_dates, k, remove_devs, ndevs)
for(i in 1:length(comp)) {
  c <- comp[[i]]
  # If enter after our selected (unshifted) waves,
  # shift to accomodate the new numbering
  if(c[1] > selected_waves[1]) {
    c[1] <- c[1] - (selected_waves[1] - 1)
  }
  # If enter before unshifted waves, set to 1
  if(c[1] < selected_waves[1]){ 
    c[1] <- 1
  }
  # If you leave after the end of unshifted waves,
  # set the end to the end of selected waves
  if(c[2] > selected_waves[length(selected_waves)]) {
    c[2] <- selected_waves[length(selected_waves)]
  }

  if(c[2] > selected_waves[1] & selected_waves[1] > 1) {
    c[2] <- c[2] - (selected_waves[1] - 1)
  }
  comp[[i]] <- c
}
changes <- sienaCompositionChange(comp)

email_cumulative_count_log_varcovar <- varCovar(log(getEmailCumulativeCountData(wave_dates, remove_devs, ndevs) + 1)[,selected_waves])
email_cumulative_count_log_cocovar <- coCovar(log(getEmailCumulativeCountData(wave_dates, remove_devs, ndevs) + 1)[,selected_waves][,1])

commit_beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, "total_commits_so_far", remove_devs)
commit_beh <- sienaDependent(commit_beh_mat[,selected_waves], type="behavior")

ownership_sm_beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, "ownership_count_simple_majority", remove_devs)
ownership_sm_beh <- sienaDependent(ownership_sm_beh_mat[,selected_waves], type="behavior")
ownership_sm_beh_varcovar <- varCovar(ownership_sm_beh_mat[,selected_waves])

ownership_count_70_beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, "ownership_count_70", remove_devs)
ownership_count_70_beh <- sienaDependent(ownership_count_70_beh_mat[,selected_waves], type="behavior")
ownership_count_70_beh_varcovar <- varCovar(ownership_count_70_beh_mat[,selected_waves])

minority_contributor_count_beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, "ownership_count_100", remove_devs)
minority_contributor_count_beh_varcovar <- varCovar(minority_contributor_count_beh_mat[,selected_waves])

ownership_count_100_beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, "ownership_count_100", remove_devs)
ownership_count_100_beh_varcovar <- varCovar(ownership_count_100_beh_mat[,selected_waves])

commit_log_mat <- log(retrieveDataMatrixFromAllData(all_data, ndevs, "total_commits_so_far", remove_devs)+1)
commit_log_varcovar <- varCovar(commit_log_mat[,selected_waves])

ownership_sm_log_mat <- log(retrieveDataMatrixFromAllData(all_data, ndevs, "ownership_count_simple_majority", remove_devs)+1)
ownership_sm_log_varcovar <- varCovar(ownership_sm_log_mat[,selected_waves])

ownership_count_70_log_mat <- log(retrieveDataMatrixFromAllData(all_data, ndevs, "ownership_count_70", remove_devs)+1)
ownership_count_70_log_varcovar <- varCovar(ownership_count_70_log_mat[,selected_waves])

ownership_count_100_log_mat <- log(retrieveDataMatrixFromAllData(all_data, ndevs, "ownership_count_100", remove_devs) + 1)
ownership_count_100_log_varcovar <- varCovar(ownership_count_100_log_mat[,selected_waves])

minority_contributor_count_log_mat <- log(retrieveDataMatrixFromAllData(all_data, ndevs, "minority_contributor_count", remove_devs) + 1)
minority_contributor_count_log_varcovar <- varCovar(minority_contributor_count_log_mat[,selected_waves])

age_mat <- retrieveAgeMatrixAllOnes(all_data, first_commit_data, ndevs, remove_devs)
age_varcovar <- varCovar(age_mat[,selected_waves])

age_raw_mat <- retrieveAgeMatrix(all_data, first_commit_data, ndevs, remove_devs)
age_log_varcovar <- varCovar(log(age_raw_mat+1)[,selected_waves])

if(beh_name != "commit_window") {
  beh_mat <- processBehNotNormalizedDataWithLeavers(all_data, first_commit_data, last_commit_data, ndevs, beh_name, remove_devs)
} else {
  beh_mat <- as.matrix(generate_commit_window_data(wave_dates, rate_den, commits, ndevs, remove_devs, "beh"))
}
beh <- sienaDependent(beh_mat[,selected_waves], type="behavior")

coevolutionData <- sienaDataCreate(changes, email_network, beh, email_cumulative_count_log_varcovar, ownership_count_70_beh_varcovar, ownership_count_100_beh_varcovar, ownership_sm_beh_varcovar, minority_contributor_count_log_varcovar, age_log_varcovar, age_varcovar, ownership_sm_log_varcovar, ownership_count_70_log_varcovar, ownership_count_100_log_varcovar, commit_log_varcovar)
coevolutionEffects <- getEffects(coevolutionData)

# Print a preliminary report of the data
print01Report(coevolutionData, coevolutionEffects, modelname="INTIAL_REPORT_MODEL_NAME_PLACEHOLDER")

# Initial effects object
# Contains example shortnames
coevolutionEffects <- includeEffects(coevolutionEffects, transTies, cycle3, outActSqrt, inPopSqrt)

# Add email and age effects
coevolutionEffects <- includeEffects(coevolutionEffects, altX, interaction1="email_cumulative_count_log_varcovar")
coevolutionEffects <- includeEffects(coevolutionEffects, altX, interaction1="age_log_varcovar")

# Add homophily
coevolutionEffects <- includeEffects(coevolutionEffects, sameX, interaction1="beh")

# Add influence
coevolutionEffects <- includeEffects(coevolutionEffects, avAlt, interaction1="email_network", name="beh")

# Create the siena algorithm object with various arguments
coevolutionAlg <- sienaAlgorithmCreate(projname="SIMULATION_PROJECT_NAME_PLACEHOLDER", diagonalize=1.0, cond=F, n3=3000, firstg=.03)

# 8 nodes in parallel, possible previous answer for guidance of the simulation
ans <- siena07ToConvergence(alg=coevolutionAlg, dat=coevolutionData, eff=coevolutionEffects, nbrNodes=8, prev=NULL)

# Print results to HTML for easier reading
siena.table(ans, type="html", file="HTML_OUTPUT_PLACEHOLDER", sig=T)

# Save our answer for future usage
save(ans, file="SIMULATION_RESULTS_PLACEHOLDER")

# Goodness of fit tests
gof.bd <- sienaGOF(ans, BehaviorDistribution, varName="beh", verbose=T, join=T, cumulative=F)
gof.id <- sienaGOF(ans, IndegreeDistribution, varName="email_network", verbose=T, join=T, levls=0:max(unlist(max_indegs)[selected_waves]))
gof.od <- sienaGOF(ans, OutdegreeDistribution, varName="email_network", verbose=T, join=T, levls=0:max(unlist(max_outdegs)[selected_waves]))
gof.gd <- sienaGOF(ans, GeodesicDistribution, varName="email_network", verbose=T, join=T)