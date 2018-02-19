get_stepwise_result <- function(pan_table, lines_order){
    # given a pan genome table and a vector with the order of lines to add, returns:
    # n_lines   n_genes
    #    1       30000
    #    3       31000
    # ...
    n_lines_vec <- c()
    n_genes_vec <- c()
    result_pav_vec <- rep(F,nrow(pan_table))
    for (n_lines in 1:length(lines_order)){
        line_to_add <- lines_order[n_lines]
        line_pav_vec <- pan_table[,line_to_add]
        result_pav_vec <- result_pav_vec | line_pav_vec
        n_genes <- sum(result_pav_vec)
        n_lines_vec <- c(n_lines_vec,n_lines)
        n_genes_vec <- c(n_genes_vec,n_genes)
    }
    return(data.frame(n_lines_vec, n_genes_vec))
}

iterate_stepwise_analysis <- function(pan_table, line_names_vec, N){
    # itearte analysis N times with randomized order
    # of line names. Returns a single df agregating results
    final_result <- data.frame(n_lines=character(), n_genes=character())
    for (x in 1:N){
        rand_lines_ord <- sample(line_names_vec)
        rand_ord_result <- get_stepwise_result(pan_table,rand_lines_ord)
        final_result <- rbind(final_result,rand_ord_result)
    }
    return(final_result)
}

# taken from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                       c(N    = length2(xx[[col]], na.rm=na.rm),
                         mean = mean   (xx[[col]], na.rm=na.rm),
                         sd   = sd     (xx[[col]], na.rm=na.rm)
                       )
                   },
                   measurevar
    )
    
    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}
