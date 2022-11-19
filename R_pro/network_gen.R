require("combinat")

# Shift the data frame up into offset, \\
# putting df[offset,] at the end of the new df
shift <- function(df, offset) {
  if (offset == nrow(df)) {
    return(df)
  }
  toshift <- data.frame(x = numeric(length(offset:nrow(df))),
  y = numeric(length(offset:nrow(df))))
  toshift[1 : (nrow(toshift) - 1), ] <- df[(offset + 1) : nrow(df), ]
  toshift[nrow(toshift), ] <- data.frame(0, 0)
  df[offset:nrow(df), ] <- toshift
  return(df)
}

# Input project name here e.g. ant
project_name <- "ant"

pro_abs_path <- "E:\\Graduate\\Trustie\\Complex Networks\\R_pro"

emails_file <- paste(pro_abs_path, "data", "ant_emails.csv", sep = "\\")
emails <- read.csv(emails_file, header = FALSE, stringsAsFactors = FALSE)
names(emails) <- c("from", "to", "date")
emails$date <- as.POSIXct(emails$date)

commits_file <- paste(pro_abs_path, "data", "ant_commits.csv", sep = "\\")
commits <- read.csv(commits_file, header = TRUE, stringsAsFactors = FALSE)
commits$committer_dt <- as.POSIXct(commits$committer_dt)

ndevs <- max(commits$committer_id, emails$from, emails$to)  # num. of devs

no_decay <- function(len, crap) {
  rep(99999, len)
}

# *****************************************
# *** Decay at some function iterations ***
# *****************************************

# Tests multiple decay rates from 1 to 100
for (rate_den in seq(1, 100, by = 1)) { # nolint
  rate <- 1 / rate_den
  decay_function <- rexp
  # Can set to no_decay for testing
  # by setting "decay_function <- no_decay"
  plots <- list()
  start_date <- min(as.POSIXct(emails$date))
  end_date <- max(as.POSIXct(emails$date))
  # Calculated this way to make sure we have every data point
  # (otherwise we may cut out some points at the end)
  ndays <- as.numeric(ceiling(end_date - start_date))
  # Change to the maximum number of iterations
  # Each iteration will be one of these units
  # e.g. if iterations_max <- ndays, each iteration
  # of the loop will be 1 day.
  iterations_max <- ndays
  # Change to the units at which time increases
  # e.g. if by day, then one_unit <- 1*60*60*24
  one_unit <- 1 * 60 * 60 * 24
  net <- matrix(0, nrow = ndevs, ncol = ndevs)
  nets <- list(matrix(0, nrow = ndevs, ncol = ndevs))
  npairs <- nrow(unique(data.frame(
                 emails$from, emails$to)))  # compute unique pairs
  # Generate pairs
  pairs <- expand.grid(seq(ndevs), seq(ndevs))
  # Remove self pairs
  pairs <- pairs[-seq(1, nrow(pairs), by = ndevs + 1), ]
  forget_time <- data.frame(pairs[, 1], pairs[2], rep(0, nrow(pairs)))
  names(forget_time) <- c("from", "to", "remaining_units")
  order <- data.frame(from = numeric(npairs), to = numeric(npairs))
  iter_count <- 0
  net_count <- 0
  prev_date <- start_date
  for (iteration in 1:iterations_max) {
    cur_date <- prev_date + one_unit
    cur_rows <- emails[which(emails$date >= prev_date
    & emails$date < cur_date), , ]
    # Update remaining units
    # forget update
    forget_time$remaining_units <- forget_time$remaining_units - one_unit
    # process emails, setting new forget_times as necessary
    if (nrow(cur_rows) > 0) {
      # set the seed here so we can test this thing
      # without random results.
      # we will remove this later when the time comes
      set.seed(iteration * rate_den)
      decay_list <- decay_function(nrow(cur_rows), rate)
      for (i in 1:nrow(cur_rows)) { # nolint
        row <- cur_rows[i, , ]
        from <- as.numeric(row$from)
        to <- as.numeric(row$to)
        date <- row$date
        # Set the forget time since an email has just been made.
        # Note that even if it's not time for the existing edge to
        # decay, this will set a brand new forget time on the edge.
        # The model is that every time an email is sent, a new forgetfulness
        # factor is attributed to the edge, regardless of the previous value.
        # Note that we subtract (cur_date - date) \\
        # because we're iterating by one_unit, and
        # the email can occur between one_unit time steps.
        # What is the meaning of "forget_time"?
        forget_time[forget_time$from == from & forget_time$to == to,
        "remaining_units"] <- decay_list[i] * one_unit - (cur_date - date)
      }
    }
    # update net every time step, i.e., one_unit
    for (i in 1:nrow(forget_time)) { # nolint
      row <- forget_time[i, ]
      net[row$from, row$to] <- 0
      if (row$remaining_units > 0) {
        net[row$from, row$to] <- 1
      }
    }
    nets[[iteration]] <- net  # nets contain every net at every time step
    prev_date <- cur_date
  }
  model_name <- paste("rate_den_", rate_den, sep = "")
  gen_data_file_path <- paste(pro_abs_path, "net_gen_data",
                              project_name, model_name, sep = "\\")
  if (!dir.exists(gen_data_file_path)) {
    dir.create(gen_data_file_path)
  }
  save(nets, file = paste(gen_data_file_path, "\\all_nets.RData", sep = ""))
  ######################################
  # The stuff for figuring out when    #
  # we lose 25%, 50%, 75% of our edges #
  # after reached the max degree.      #
  ######################################
  decay_ind <- function(dev, percent, max_inds, mat, max_degs) {
    which_inds <- which(mat[, dev] <= floor(max_degs[dev] * percent))
    which_inds[which(which_inds > max_inds[dev])][1]
  }
  outdegs <- lapply(nets, rowSums)
  indegs <- lapply(nets, colSums)
  # Each row is a time frame, each column is a dev
  outdegs_mat <- do.call(rbind, lapply(outdegs, unlist))
  indegs_mat <- do.call(rbind, lapply(indegs, unlist))
  max_outdegs_inds <- apply(outdegs_mat, 2, which.max)
  max_indegs_inds <- apply(indegs_mat, 2, which.max)
  max_outdegs <- apply(outdegs_mat, 2, function(x) x[which.max(x)])
  max_indegs <- apply(indegs_mat, 2, function(x) x[which.max(x)])
  # The below are descriptive statistics that may
  # be interesting to look at.
  # They're can be commented out if you're
  # not interested in descriptive statistics.
  print("doing degs")
  # outdegs_#_inds and indegs_#_inds
  outdegs_75_inds <- sapply(1:ndevs, decay_ind, percent = .75,
  max_inds = max_outdegs_inds, mat = outdegs_mat, max_degs = max_outdegs)
  indegs_75_inds <- sapply(1:ndevs, decay_ind, percent = .75,
  max_inds = max_indegs_inds, mat = indegs_mat, max_degs = max_indegs)
  outdegs_50_inds <- sapply(1:ndevs, decay_ind, percent = .50,
  max_inds = max_outdegs_inds, mat = outdegs_mat, max_degs = max_outdegs)
  indegs_50_inds <- sapply(1:ndevs, decay_ind, percent = .50,
  max_inds = max_indegs_inds, mat = indegs_mat, max_degs = max_indegs)
  outdegs_25_inds <- sapply(1:ndevs, decay_ind, percent = .25,
  max_inds = max_outdegs_inds, mat = outdegs_mat, max_degs = max_outdegs)
  indegs_25_inds <- sapply(1:ndevs, decay_ind, percent = .25,
  max_inds = max_indegs_inds, mat = indegs_mat, max_degs = max_indegs)
  # diff_outdegs and diff_indegs (i.e., negative value of days
  # when we loss 25%, 50% and 75% edges)
  diff_outdegs_75 <- max_outdegs_inds - outdegs_75_inds
  diff_outdegs_50 <- max_outdegs_inds - outdegs_50_inds
  diff_outdegs_25 <- max_outdegs_inds - outdegs_25_inds
  diff_indegs_75 <- max_indegs_inds - indegs_75_inds
  diff_indegs_50 <- max_indegs_inds - indegs_50_inds
  diff_indegs_25 <- max_indegs_inds - indegs_25_inds
  # save datas
  save(diff_outdegs_75, file = paste(gen_data_file_path,
       "\\diff_outdegs_75.RData", sep = ""))
  save(diff_outdegs_50, file = paste(gen_data_file_path,
       "\\diff_outdegs_50.RData", sep = ""))
  save(diff_outdegs_25, file = paste(gen_data_file_path,
       "\\diff_outdegs_25.RData", sep = ""))
  save(diff_indegs_75, file = paste(gen_data_file_path,
       "\\diff_indegs_75.RData", sep = ""))
  save(diff_indegs_50, file = paste(gen_data_file_path,
       "\\diff_indegs_50.RData", sep = ""))
  save(diff_indegs_25, file = paste(gen_data_file_path,
       "\\diff_indegs_25.RData", sep = ""))
  # Can add more into here for different wave split counts
  # Split emails to k parts.
  for (k in c(8)) {
    split_ind <- ceiling(seq(nrow(emails) / k,
    nrow(emails), by = nrow(emails) / k))
    split_dates <- emails[split_ind, "date"]
    inds <- c()
    for (date in split_dates) {
      inds <- c(inds, which(abs(as.numeric(
        seq(start_date, end_date, by = one_unit)) - date)
        == min(abs(as.numeric(
        seq(start_date, end_date, by = one_unit)) - date))))
    }
    # model_name <- paste(project_name, "_", rate_den, sep = "") # nolint
    # The sum of emails between two savenets is equal 
    # and equal to nrow(emails) / k.
    savenets <- nets[inds]
    save(savenets, file = paste(gen_data_file_path, "\\",
                                k, "split_nets.RData", sep = ""))
  }
}
