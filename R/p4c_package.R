
#' Create a new p4cProfile object.
#' 
#' \code{p4cNewProfile} constructs a p4cProfile object from a genomic 2D track.
#' 
#' This function is the primary object constructor. It processes a 2D genomic track 
#' with (UMI)-4C raw contact data and returns a p4cProfile object. 
#' 
#' @param track_nm 2D genomic track name. Can also be a character vector of track
#'  names that will be pooled together by summation. 
#' @param bait_chrom Character. Bait chromosome.
#' @param bait_start Integer. Start coordinate of the bait.
#' @param scope_5,scope_3 Define the scope of the profile. How many bases upstream
#'  and downstream of the bait should be used for generating the profile. 
#'  Default is set to 250kb upstream and downstream.
#' @param stat_type Define the method for calculating the profile contact intensities. 
#'  Either "linear" or "log", for arithmetic or geometric means, respectively. Default 
#'  is "linear".
#' @return Returns a p4cProfile object.
#' 
#' @export
p4cNewProfile <- function(track_nm, bait_chrom = gtrack.attr.get(track_nm[1], "Bait_chr"), 
    bait_start = as.numeric(gtrack.attr.get(track_nm[1], "Bait_coord")), scope_5=2.5e5, 
    scope_3=2.5e5, stat_type)
    {
    for (nm in track_nm)
    {
        if (!gtrack.exists(nm))
        {
            stop("cannot find track ", nm)
        }
    }
    if (nchar(bait_chrom) < 1 | nchar(bait_start) < 1 | is.na(bait_start))
    {
        stop("No bait attributes for track name - Please supply the parameters")
    }
    
    .checkConf()
    
    # Check that the scopes are compatible with bait and chromosome size
    if (bait_start < scope_5){
        warning("The bait coordinate is", bait_start, "but the scope_5 was set to", scope_5,
                ".\nThis would give negative coordinate. scope_5 is adjusted.")
        scope_5 <- bait_start
    }
    chrom_size <- gintervals.all()[gintervals.all()$chrom == bait_chrom, ]$end
    if (bait_start + scope_3 > chrom_size){
        warning("The bait coordinate is", bait_start, "and the chromosome size is",
                chrom_size, "but the scope_3 was set to", scope_3,
                ".\nThis would give coordinate above chromsome size. scope_3 is adjusted.")
        scope_3 <- chrom_size - bait_start
    }    
    p4c_obj <- list(track_nm = paste(track_nm, collapse = " "), bait = list(chrom = bait_chrom, 
        start = as.numeric(bait_start), bait_lookup_expansion = as.numeric(getOption("TG3C.bait_lookup_expansion"))), 
        scope = c(scope_5, scope_3))
    
    class(p4c_obj) <- "p4cProfile"
    names(p4c_obj$scope) <- c("scope_5", "scope_3")
    if (missing(stat_type))
    {
        p4c_obj <- .p4cSetDgramParams(p4c_obj)
    } else
    {
        p4c_obj <- .p4cSetDgramParams(p4c_obj, stat_type = stat_type)
    }
    
    p4c_obj <- .p4cSetDgramGraphics(p4c_obj)
    p4c_obj <- .p4cSetPngParams(p4c_obj)
    
    stat_type <- p4c_obj$dgram_params$stat_type
    p4c_obj <- .p4cGenerateBaitDgram(p4c_obj)
    return(p4c_obj)
}


#### Generic function for print

#' @export
print.p4cProfile <- function(p4c_obj)
{
    cat("p4cProfile object\n\n")
    cat("4C Track name: ", p4c_obj$track_nm, "\n\n")
    cat("Current bait: \n")
    cat("Chrom:", p4c_obj$bait$chrom, " Start:", p4c_obj$bait$start, "\n\n")
    cat("Current scope: \n")
    cat("scope_5:", p4c_obj$scope["scope_5"], " scope_3:", p4c_obj$scope["scope_3"], 
        "\n")
    cat("stat type:", p4c_obj$dgram_params$stat_type, "\n")
}


#### Summery generic function will return basic stats of track
#' @export
summary.p4cProfile <- function(p4c_obj)
{
    track_nm <- unlist(strsplit(p4c_obj$track_nm, " "))[1]
    bait_chrom <- paste0('chr', gsub('chr', '', p4c_obj$bait$chrom))
    bait_start <- p4c_obj$bait$start
    bait_lookup_expansion <- p4c_obj$bait$bait_lookup_expansion
    
   
    d <- nrow(ALLGENOME[[1]])
    ivals <- gintervals.2d(rep(bait_chrom, d), rep(bait_start - bait_lookup_expansion, d), rep(bait_start + 
        bait_lookup_expansion, d), ALLGENOME[[1]]$chrom, ALLGENOME[[1]]$start, ALLGENOME[[1]]$end)
    
    ivals.trans <- ivals[ivals$chrom2 != bait_chrom, ]
    
    ivals.nearcis <- gintervals.2d(bait_chrom, bait_start - bait_lookup_expansion, bait_start + 
        bait_lookup_expansion, bait_chrom, bait_start - 1e+06, bait_start + 1e+06)
    filter_arg <- paste0(track_nm, " > 0")
    track_summary <- gsummary(track_nm, intervals = gscreen(filter_arg, ivals))
    track_trans_summary <- gsummary(track_nm, intervals = gscreen(filter_arg, ivals.trans))
    track_nearcis_summary <- gsummary(track_nm, intervals = gscreen(filter_arg, ivals.nearcis))
    
    covered_fends = unname(track_summary["Total intervals"])
    total_contacts = unname(track_summary["Sum"])
    cis_contacts = total_contacts - unname(track_trans_summary["Sum"])
    cis_ratio = cis_contacts/total_contacts
    nearcis_contacts = unname(track_nearcis_summary["Sum"])
    range = track_summary[c("Min", "Max")]
    
    vec <- c(`Covered fends` = round(covered_fends), `Total contacts` = round(total_contacts), 
        `Cis contacts` = round(cis_contacts), `Cis ratio` = round(cis_ratio, 2), 
        `Nearcis (<1Mb) contacts` = round(nearcis_contacts))
    vec
}

#' Plot 4C nearcis profile.
#' 
#' Plot a nearcis profile in the defined genomic scope. Can be applied on a single 
#' p4cProfile object or a two object for a comparison plot.
#' 
#' This function calls two methods depending on whether it is called on a single object
#' or two objects. When it is called on a single p4cProfile object it will produce a nearcis
#' plot with smoothed trendline, and a contact intensity domainogram. 
#' Calling this function for two object will produce a comparison plot of where the 
#' 4C profiles are normalized for UMI coverage. The two p4cProfile objects must be
#' constructed for the \strong{same bait} with the \strong{same genomic scope}.
#' 
#' @param p4c_obj1 p4cProfile object created by \code{p4cNewProfile}.
#' @param p4c_obj2 Optional. A second p4cProfile object that will be compared to p4c_obj1.
#' @param trend_scale This parameter controls the smoothing method for the trendline. Either 
#'  'adaptive' (default), or, for fixed window trend, an integer with the desired size 
#'  of the window (in restriction fragments units).
#' @param min_win_cov If the default smoothing method is used ('adaptive'), this parameter
#'  controls the window size by requiring that no less than \code{min_win_cov} molecules will be included
#'  in that genomic window.
#' @param png_fn File to save the plot to. The default (NULL) will result in plotting on the 
#'  current graphic device.
#' @param plot.colorbar Optional. Logical that defines whether a colorbar for the domainogram 
#'  should be plotted.
#' @param xlim Similar to R basic graphics. A two elements vector that controls the genomic 
#'  focus of the plot (X axis). 
#' @param ylim Similar to R basic graphics. A two elements vector that controls the trendline
#'  ylim (In units of molecules).
#' @param main Plot title.
#' @param sd Controls the margins of the trendline. The Number of bionomial theoretical standard deviations.
#' @param legend.text Character vector with length 2. Controls the text in the color legend of the trendline.
#'  The default are the names of the profiles.    
#' @param col Vector with length 2 that controls the colors of each profile in a comparative plot. 
#' @param trend_only Optional. Logical that defines whether only the trend plot should be displayed.
#' @param filename1 Optional. File for the bedgraph with the trendline of `p4c_obj1`
#' @param filename2 Optional. File for the bedgraph with the trendline of `p4c_obj2`
#' 
#' @examples 
#' \donttest{
#' # Create a p4cProfile object:
#' fc <- p4cNewProfile("umi4C_example_CMK_ANK1_TSS")
#' 
#' # Plot a profile on:
#' plot(fc, xlim=c(41554693, 41754693), ylim=c(0, 5), plot.colorbar=TRUE)
#' 
#' 
#' # Comparing two profiles:
#' # Create a second profile
#' fc_cond2 <- p4cNewProfile("umi4C_example_DND41_ANK1_TSS")
#' 
#' # Plot profiles
#' plot(fc, fc_cond2, col=c('red', 'blue'), min_win_cov=100, main = '4C figure', png_fn='fig.png') 
#' }
#' 
#' @export
plot.p4cProfile <- function(p4c_obj1, p4c_obj2, trend_scale = "adaptive", png_fn = NULL, 
    ...)
    {
    if (missing(p4c_obj2))
    {
        plotSingleProf(p4c_obj = p4c_obj1, trend_scale = trend_scale, png_fn = png_fn, 
            ...)
    }
    
    if (!missing(p4c_obj1) & !missing(p4c_obj2))
    {
        plotCompProf(p4c_obj1, p4c_obj2, trend_scale = trend_scale, png_fn = png_fn, 
            ...)
    }
    
}

#' @rdname plot.p4cProfile
#' @export
plotSingleProf <- function(p4c_obj, trend_scale, png_fn, plot.colorbar, add.func, 
    xlim, ylim, trend_only, main, sd, ...) UseMethod("plotSingleProf")

#' @export
plotSingleProf.p4cProfile <- function(p4c_obj, png_fn = NULL, trend_scale = "adaptive", 
    ref_track_nm = NA, min_win_cov = 50, plot.colorbar = FALSE, add.func, xlim, ylim, 
    trend_only = FALSE, main, sd = 2, filename1, ...)
    {
    if (trend_scale != "adaptive" & !is.numeric(trend_scale))
    {
        stop("illegal trend scale ", trend_scale)
    }
    
    if (is.numeric(trend_scale) & !trend_scale %in% p4c_obj$dgram_params$dgram_scales)
    {
        stop("trend_scale was not generated for ", trend_scale, ": set dgram_scales in p4c.conf")
    }
    
    if (is.numeric(trend_scale))
    {
        # the true index of the desired trend_scale in the dgram matrix
        trend_scale_idx <- which(p4c_obj$dgram_params$dgram_scales == trend_scale) + 
            2
    }
    
    if (is.null(p4c_obj[["dgram"]]))
    {
        warning("No dgram matrix found for p4C obj, please create one using .p4cGenerateBaitDgram()\n
            created a temporary matrix")
        p4c_obj <- .p4cGenerateBaitDgram(p4c_obj)
    }
    
    dgram <- p4c_obj$dgram  #access dgram table
    coords <- dgram[, 1]
    n <- length(dgram[, 1])
    
    bait_x <- p4c_obj$bait$start
    if (missing(xlim))
    {
        horiz5 <- bait_x - p4c_obj$scope[1]
        horiz3 <- bait_x + p4c_obj$scope[2]
    } else
    {
        horiz5 <- xlim[1]
        horiz3 <- xlim[2]
    }
    
    if (is.numeric(trend_scale))
    {
        coords <- .p4cWinGeoMeanCoordinate(coords, trend_scale * 2, bait_x)
        coords <- coords[trend_scale:(n - trend_scale)]
        dgram <- dgram[trend_scale:(n - trend_scale), ]  #truncate the margins of the dataframe to remove NAs
        trend <- dgram[, trend_scale_idx]
        cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (fixed trend scale: %i fragments)",
            cur_name, trend_scale)
        track_name <- paste0(cur_name,"_smoothed",trend_scale,"frags")
        message("fixed trend scale: ", trend_scale)
    } else if (trend_scale == "adaptive")
    {
        if ("smoothedTrend" %in% names(p4c_obj))
        {
            message("trend is adaptive")
        } else
        {
            message("Will create smoothed trend for track with min_win_cov = ", min_win_cov)
            p4c_obj <- .p4cSmoothedTrend(p4c_obj, min_win_cov)
        }
        coords <- p4c_obj$smoothedTrend$start
        trend <- p4c_obj$smoothedTrend$trend
        cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (min coverage - %i molecules)", 
            cur_name, min_win_cov)
        track_name <- paste0(cur_name,"_smoothedMin",min_win_cov,"mols")
    }
    
    trend <- trend[!is.na(coords)]
    coords <- na.omit(coords)

    if(!missing(filename1)){
        cur_name <- p4c_obj$track_nm
        cur_col <- paste((col2rgb(col)), collapse = ",")
        header_line1 <- sprintf("browser position %s:%i-%i", paste0("chr", gsub("chr", 
            "", p4c_obj$bait$chrom)), coords[1], tail(coords, 1))
        header_line2 <- sprintf("track type=bedGraph name=\"%s\" description=\"%s\" color=%s", 
            track_name, cur_desc, cur_col)
        write(header_line1, filename1)
        write(header_line2, filename1, append = T)
    
        bed_df <- data.frame(chrom = paste0("chr", gsub("chr", "", p4c_obj$bait$chrom)), 
            start = round(coords), end = round(coords) + 1, 
            value = trend)
    
        write.table(bed_df, filename1, append = T, quote = F, col.names = F, row.names = F, sep="\t")
        message("Wrote bedGraph file for ",cur_name," to ", filename1)
    }
    
    if (!is.null(png_fn))
    {
        png_fn <- paste0(p4c_obj$figure_params$png_dir, "", png_fn)
        png(png_fn, width = p4c_obj$figure_params$png_w, height = p4c_obj$figure_params$png_h, 
            res = p4c_obj$figure_params$png_res)
    }
    
    par(mar = c(0, 4, 1.2, 2))
    par(oma = c(2, 0.4, 1, 0.4))
    par(xaxs = "i", yaxs = "i")
    
    if (plot.colorbar)
    {
        layout(matrix(c(1, 2, 3, 3), 2, 2), heights = c(3.2, 2.8), c(12, 1))
    }
    
    if (!trend_only & !plot.colorbar)
    {
        layout(matrix(c(1, 2), 2, 1), heights = c(3.2, 2.8))
    }
    
    # remove from plot long runs (>30kb) of missing RE site
    coords_fends <- dgram[, 1]
    missing_re_idx <- which(coords_fends[-1] - coords_fends[-length(coords_fends)] > 
        30000)
    message("Plotting from ", horiz5, " to ", horiz3)
    
    if (missing(ylim))
    { 
        ylim <- c(0, (0.8 + max(trend, na.rm = T)))
    }
    
    # plot trend
    plot(coords, pmin(trend, ylim[2]), type = "l", lwd = 2, col = "black", axes = F, 
        ylab = "", ylim = ylim, xlim = c(horiz5, horiz3))
    
    # add title
    if (!missing(main)) 
    {
        title(main = main)
    } 
    
    # add y axis
    axis(2, las = 1, cex.axis = 1)
    
    omitted_idx <- attr(coords, "na.action")
    
    # plot discrete data points
    stat_type <- p4c_obj$dgram_params$stat_type
    if (stat_type == "log") 
        points(dgram[, 1], log2(1 + dgram[, 2]), col = "gray", pch = 19, cex = 0.5)
    if (stat_type == "linear") 
        points(dgram[, 1], pmin(dgram[, 2], ylim[2]), col = "gray", pch = 19, cex = 0.5)
    
    # plot sd of trend sqrt(k/n)
    if (is.numeric(trend_scale))
    {
        trend_nna <- trend
        trend_nna[is.na(trend_nna)] <- 0
        message("trend_scale = ", trend_scale)
        polygon(c(coords, rev(coords)), c(trend_nna - sd * sqrt(trend_nna/(trend_scale * 
            2)), rev(trend_nna + sd * sqrt(trend_nna/(trend_scale * 2)))), col = "grey", 
            border = 0)
    }
    if (trend_scale == "adaptive")
    {
        trend_nna <- trend
        trend_nna[is.na(trend_nna)] <- 0
        scales <- p4c_obj$smoothedTrend$scale
        scales[is.na(scales)] <- 1
        polygon(c(coords, rev(coords)), c(trend_nna - sd * sqrt(trend_nna/(scales * 
            2)), rev(trend_nna + sd * sqrt(trend_nna/(scales * 2)))), col = "grey", 
            border = 0)
    }
    
    # plot main trend
    lines(coords, trend, type = "l", lwd = 3, col = "black")
    
    if (trend_only){
        if (!is.null(png_fn))
        {
            dev.off()
        }
        return()
    }
    
    # plot dgram
    par(mar = c(4, 4, 0, 2))
    par(mgp = c(2.5, 0.8, 0))
    par(xaxs = "i", yaxs = "i")
    
    n = ncol(dgram)
    dgram_orig = dgram
    norm_factor = max(dgram[,3:n], na.rm = T)
    dgram[, 3:n] = log2(dgram[, 3:n]) - log2(norm_factor)
    if (norm_factor < 1)
    {
        dgram[, 3:n] = log2(dgram[, 3:n] + 1) - log2(norm_factor + 1)  #fixed a bug that resulted in strange output if max value of trend <1
    }
    message("Range of values in dgram: from ", min(dgram_orig[, n:3], na.rm = T), 
        " to ", max(dgram_orig[, n:3], na.rm = T))
    message("Transformed to log2(fraction of max): from ", min(dgram[, n:3], na.rm = T), 
        " to ", max(dgram[, n:3], na.rm = T))
    shades <- colorRampPalette(p4c_obj$graphics_params$shades)(1000)
    
    image(x = coords, z = pmin(dgram[, n:3], p4c_obj$graphics_params$zlim[2]), col = shades, 
        xlab = sprintf("chrom %s coordinates (Mb)", gsub("chr", "", p4c_obj$bait$chrom)), 
        zlim = p4c_obj$graphics_params$zlim, xlim = c(horiz5, horiz3), cex.lab = 1.2, 
        yaxt = "n", xaxt = "n")
    
    axis(1, dgram[, 1], labels = F, col.ticks = "light grey")
    axis(1, round(seq(horiz5, horiz3, l = 9)/1e+06, 2) * 1e+06, labels = round(seq(horiz5, 
        horiz3, l = 9)/1e+06, 2), cex.lab = 2, cex.axis = 1.4)
    
    scales_l <- length(p4c_obj$dgram_params$dgram_scales)
    axis(2, c(0, 0.5), labels = c(p4c_obj$dgram_params$dgram_scales[scales_l], 
        p4c_obj$dgram_params$dgram_scales[scales_l/2]), 
        las = 1, cex.axis = 1)
    
    # mask missing re in dgram
    sapply(missing_re_idx, function(i) polygon(c(coords_fends[i] - 10000, coords_fends[i], 
        coords_fends[i + 1], coords_fends[i + 1] + 10000), c(0, 1, 1, 0), col = "white", 
        border = NA))
    
    # plot colorbar
    if (plot.colorbar)
    {
        par(mar = c(8, 4, 8, 4))
        par(mgp = c(2, 0.8, 0))
        .image.scale(zlim = p4c_obj$graphics_params$zlim, col = colorRampPalette(p4c_obj$graphics_params$shades)(100), 
            horiz = F, axes = F, xlab = "", ylab = "", las = 1)
        axis_range <- p4c_obj$graphics_params$zlim
        expand_axis_range <- seq(axis_range[1], axis_range[2], length = 5)
        axis(2, at = expand_axis_range, labels = sprintf("%.02f%%", 100 * (2^expand_axis_range)), 
            las = 2)
        title(sub = "contacts/window")
    }
    # execute a function after plotting
    if (!missing(add.func))
    {
        eval(add.func)
    }
    if (!is.null(png_fn))
    {
        dev.off()
    }
    
}

#' @rdname plot.p4cProfile
#' @export
plotCompProf <- function(p4c_obj1, ref_p4c_obj, trend_scale, png_fn, col, min_win_cov, 
    xlim, zlim, legend.text, ylim, dgram.method, main, sd) UseMethod("plotCompProf")

#' @export
plotCompProf.p4cProfile <- function(p4c_obj1, ref_p4c_obj, trend_scale = "adaptive",
        png_fn = NULL, col = c("red", "black"), 
        min_win_cov = 50, xlim, zlim = c(-1.5, 1.5), legend.text, ylim, 
        dgram.method = "delta", main, sd = 2, filename1, filename2,
        trend_only = FALSE, ...)
        {
    if (trend_scale != "adaptive"){
        message("Only the trend scale adaptive is currently supported in plot",
            " when 2 p4cProfiles are provided.")
    }
    # normalize p4c_obj1 to p4c_obj2
    p4c_obj1 <- .p4cNormDgram(p4c_obj1, ref_p4c_obj)
    
    # call adaptive smoothing function to create comparative list
    comp_list <- .p4cSmoothedTrendComp(p4c_obj1, ref_p4c_obj, min_win_cov = min_win_cov)
    
    coords <- comp_list$trend_mat$start
    trend1 <- comp_list$trend_mat$p4c_obj1
    ref_trend <- comp_list$trend_mat$ref_p4c_obj
    
    if(!missing(filename1)){
        cur_name <- p4c_obj1$track_nm
        cur_col <- paste((col2rgb(col[1])), collapse = ",")
        cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (min coverage - %i molecules) when norm to %s",
            cur_name, min_win_cov,ref_p4c_obj$track_nm)
        track_name <- paste0(cur_name,"_smoothedMin",min_win_cov,"mols_normTo",ref_p4c_obj$track_nm)
        header_line1 <- sprintf("browser position %s:%i-%i", paste0("chr", gsub("chr", 
            "", p4c_obj1$bait$chrom)), coords[1], tail(coords, 1))
        header_line2 <- sprintf("track type=bedGraph name=\"%s\" description=\"%s\" color=%s", 
            track_name, cur_desc, cur_col)
        write(header_line1, filename1)
        write(header_line2, filename1, append = T)
    
        bed_df <- data.frame(chrom = paste0("chr", gsub("chr", "", p4c_obj1$bait$chrom)), 
            start = round(coords), end = round(coords) + 1, 
            value = trend1)
    
        write.table(bed_df, filename1, append = T, quote = F, col.names = F, row.names = F, sep="\t")
        message("Wrote bedGraph file for ",cur_name," to ", filename1)
    }
    if(!missing(filename2)){
        cur_name <- ref_p4c_obj$track_nm
        cur_col <- paste((col2rgb(col[2])), collapse = ",")
        cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (min coverage - %i molecules)", 
            cur_name, min_win_cov)
        track_name <- paste0(cur_name,"_smoothedMin",min_win_cov,"mols_toNorm",p4c_obj1$track_nm)
        header_line1 <- sprintf("browser position %s:%i-%i", paste0("chr", gsub("chr", 
            "", p4c_obj1$bait$chrom)), coords[1], tail(coords, 1))
        header_line2 <- sprintf("track type=bedGraph name=\"%s\" description=\"%s\" color=%s", 
            track_name, cur_desc, cur_col)
        write(header_line1, filename2)
        write(header_line2, filename2, append = T)
    
        bed_df <- data.frame(chrom = paste0("chr", gsub("chr", "", p4c_obj1$bait$chrom)), 
            start = round(coords), end = round(coords) + 1, 
            value = ref_trend)
    
        write.table(bed_df, filename2, append = T, quote = F, col.names = F, row.names = F, sep="\t")
        message("Wrote bedGraph file for ",cur_name," to ", filename2)
    }
    
    if (missing(xlim)) 
        xlim = range(coords, na.rm = T)
    
    dlt_dgram <- comp_list$dlt_dgram
    n <- ncol(dlt_dgram)
    
    # plot comp profiles
    if (!is.null(png_fn))
    {
        png_fn <- paste0(p4c_obj1$figure_params$png_dir, "", png_fn)
        png(png_fn, width = p4c_obj1$figure_params$png_w, height = p4c_obj1$figure_params$png_h, 
            res = p4c_obj1$figure_params$png_res)
    }
    
    if(!trend_only)
    {
        layout(matrix(1:3, nrow = 3), heights = c(2.5, 0.2, 2))
    }
    par(cex = 1)
    par(mar = c(0, 4, 0, 2))
    par(oma = c(2, 0.4, 1, 0.4))
    par(xaxs = "i", yaxs = "i")
    
    if (missing(ylim)) 
        ylim <- c(0, max(trend1, ref_trend, na.rm = T) + 0.5)
    
    message("Plotting ", xlim[1], " ", xlim[2])
    message("Range ", range(coords)[1], " to ", range(coords)[2])
    
    plot(coords, trend1, lwd = 2, type = "l", ylim = ylim, xlab = "", ylab = "", 
        xaxt = "n", yaxt = 'n', col = col[1], xlim = xlim, bty = "n")
    
    axis(2, las = 1, cex.axis = 1)
    points(p4c_obj1$dgram.norm$dgram[, 1], pmin(p4c_obj1$dgram.norm$dgram[, 2], ylim[2]), 
        col = col[1], pch = 16, cex = 0.7)
    points(ref_p4c_obj$dgram[, 1], pmin(ref_p4c_obj$dgram[, 2], ylim[2]), col = col[2], 
        pch = 16, cex = 0.7)
    
    # sd for 1 (2 sd from the mean):
    trend_nna <- trend1
    trend_nna[is.na(trend_nna)] <- 0
    scales <- comp_list$trend_mat$scale
    scales[is.na(scales)] <- 1
    polygon(c(coords, rev(coords)), c(trend_nna - sd * sqrt(trend_nna/(scales * 2)), 
        rev(trend_nna + sd * sqrt(trend_nna/(scales * 2)))), col = "grey", border = "dark grey", 
        density = 55)
    lines(coords, trend1, lwd = 2, type = "l", col = col[1])
    
    # sd for 2:
    trend_nna <- ref_trend
    trend_nna[is.na(trend_nna)] <- 0
    scales <- comp_list$trend_mat$scale
    scales[is.na(scales)] <- 1
    polygon(c(coords, rev(coords)), c(trend_nna - sd * sqrt(trend_nna/(scales * 2)), 
        rev(trend_nna + sd * sqrt(trend_nna/(scales * 2)))), col = "grey", border = "dark grey", 
        density = 80)
    lines(coords, trend1, lwd = 5, type = "l", col = col[1])
    lines(coords, ref_trend, lwd = 5, type = "l", col = col[2])
    
    leg_text <- paste0(c(paste0(p4c_obj1$track_nm, " N=", round(sum(p4c_obj1$dgram[, 
        2], na.rm = T))), paste0(ref_p4c_obj$track_nm, " N=", round(sum(ref_p4c_obj$dgram[, 
        2], na.rm = T)))))
    
    if (!missing(legend.text))
    {
        leg_text <- legend.text
    }
    
    legend("topleft", legend = leg_text, fill = col)
    
    if (!missing(main)) 
        title(main = main)
    
    if (trend_only){
        if (!is.null(png_fn))
        {
            dev.off()
        }
        return()
    }
    
    par(mar = c(0, 4, 0, 2))
    
    # plot heatmap of smoothing
    scales <- comp_list$trend_mat$scale
    s_cols <- c("orange", "dark blue", "black")
    image(x = coords, z = as.matrix(scales), axes = F, col = colorRampPalette(s_cols)(100), 
        xlim = xlim)
    
    # plot delta dgram
    par(mar = c(0, 4, 0, 2))
    if (dgram.method == "delta")
    {
        message("Range of dlt_dgram: ", min(dlt_dgram[, n:2]), " ", max(dlt_dgram[, 
            n:2]))
        image(x = coords, z = dlt_dgram[, n:2], axes = F, col = colorRampPalette(c(col[2], 
            "white", col[1]))(1000), zlim = zlim, bty = "o", xlim = xlim, )
    } else if (dgram.method == "p4c_obj1" | dgram.method == "ref_p4c_obj")
    {
        if (dgram.method == "p4c_obj1") 
            dgram <- p4c_obj1$dgram
        if (dgram.method == "p4c_obj2") 
            dgram <- ref_p4c_obj$dgram
        
        n = ncol(dgram)
        dgram_orig = dgram
        norm_factor = max(dgram[, 3:n], na.rm = T)
        
        if (norm_factor < 1)
        {
            dgram[, 3:n] = log2(dgram[, 3:n] + 1) - log2(norm_factor + 1) 
        } else
        {
            dgram[, 3:n] = log2(dgram[, 3:n]) - log2(norm_factor)
        }
        message("Range of values in dgram: from ", min(dgram_orig[, n:3], na.rm = T), 
            " to ", max(dgram_orig[, n:3], na.rm = T))
        message("Transformed to log2(fraction of max): from ", min(dgram[, n:3], na.rm = T), 
            " to ", max(dgram[, n:3], na.rm = T))
        shades <- colorRampPalette(p4c_obj1$graphics_params$shades)(1000)
        
        image(x = coords, z = pmin(dgram[, n:3], p4c_obj1$graphics_params$zlim[2]), 
            col = shades, xlab = sprintf("chrom %s coordinate (Mb)", gsub("chr", 
                "", p4c_obj1$bait$chrom)), zlim = p4c_obj1$graphics_params$zlim, 
            xlim = xlim, cex.lab = 1.2, yaxt = "n", xaxt = "n")
    }
    
    axis(1, round(seq(xlim[1], xlim[2], l = 9)/1e+06, 2) * 1e+06, labels = round(seq(xlim[1], 
        xlim[2], l = 9)/1e+06, 2), cex.lab = 2, cex.axis = 1.4)
    if (!is.null(png_fn))
    {
        dev.off()
    }
}

# function to return data on interval
.p4cReportIntervalData <- function(p4c_obj, ref_p4c_obj, start, end, min_win_cov = 30, 
    rawdata = FALSE)
    {
    # normalize p4c_obj1 to p4c_obj2
    p4c_obj <- .p4cNormDgram(p4c_obj, ref_p4c_obj)
    
    # call adaptive smoothing function to create comparative list
    comp_list <- .p4cSmoothedTrendComp(p4c_obj, ref_p4c_obj, min_win_cov = min_win_cov)
    
    if (!rawdata)
    {
        scope_idx <- which(comp_list$trend_mat$start >= start & comp_list$trend_mat$start <= 
            end)
        
        
        slice_df <- comp_list$trend_mat[scope_idx, ]
    } else
    {
        scope_idx <- which(ref_p4c_obj$dgram[, 1] >= start & ref_p4c_obj$dgram[, 
            1] <= end)
        
        slice_df <- data.frame(p4c_obj$dgram.norm$dgram[scope_idx, 1:2], ref_p4c_obj$dgram[scope_idx, 
            2], rep(NA, length(scope_idx)))
    }
    
    colnames(slice_df) <- c("start", p4c_obj$track_nm, ref_p4c_obj$track_nm, "scale")
    slice_df$log2_foldChange <- log2(slice_df[, 2]/slice_df[, 3])
    
    # transform to fend space the scales
    slice_df$scale <- slice_df$scale * 2
    
    
    
    slice_df
}

#' Compare contacts intensities between two profiles
#' 
#' Compare contact intensities at specific genomic intervals. 
#' 
#' When a genomic region is suspected to be defferentialy contacting between two profiles, 
#' it is possible to compare these two profiles in the same intervals by deriving the normalized contact values
#' in that region. This function facilitates this task. It returns a data.frame with 
#' the either the row molecule counts if \code{rowdata = TRUE} or normalized other wise.
#' The smoothing is controlled by \code{min_win_cov}.
#' 
#' The p-value is calculated for the entire interval with a chi-square test. However,
#' it assumes the interval is <1Mb from the bait. If the interval is >1Mb from the bait
#'  \code{max_scope} should be increased accordingly.
#' 
#' @param p4c_obj A \code{p4cProfile} object.
#' @param ref_p4c_obj The \code{p4cProfile} object one wishes to compare to.
#' @param start,end Start and end coordinates of the genomic interval.
#' @param min_win_cov Controls the smoothing as described in \code{plot.p4cProfile}.
#' @param rowdata Return raw molecule counts instead of smoothed data.
#' 
#' @export
p4cIntervalsMean <- function(p4c_obj, ref_p4c_obj, start, end, min_win_cov = 30, 
    rawdata = FALSE, max_scope=1e6)
    {
    temp1 <- p4cNewProfile(unlist(strsplit(p4c_obj$track_nm, " ")), p4c_obj$bait$chrom, 
        p4c_obj$bait$start, max_scope, max_scope)
    temp2 <- p4cNewProfile(unlist(strsplit(ref_p4c_obj$track_nm, " ")), ref_p4c_obj$bait$chrom, 
        ref_p4c_obj$bait$start, max_scope, max_scope)
    
    # normalize p4c_obj1 to p4c_obj2
    p4c_obj <- .p4cNormDgram(p4c_obj, ref_p4c_obj)
    
    # call adaptive smoothing function to create comparative list
    comp_list <- .p4cSmoothedTrendComp(p4c_obj, ref_p4c_obj, min_win_cov = min_win_cov)
    
    if (!rawdata)
    {
        scope_idx <- which(comp_list$trend_mat$start >= start & comp_list$trend_mat$start <= 
            end)
        
        
        slice_df <- data.frame(min(comp_list$trend_mat[scope_idx, 1]), max(comp_list$trend_mat[scope_idx, 
            1]), mean(comp_list$trend_mat[scope_idx, 2]), mean(comp_list$trend_mat[scope_idx, 
            3]), mean(comp_list$trend_mat[scope_idx, 4]))
    } else
    {
        scope_idx <- which(ref_p4c_obj$dgram[, 1] >= start & ref_p4c_obj$dgram[, 
            1] <= end)
        
        slice_df <- data.frame(min(p4c_obj$dgram.norm$dgram[scope_idx, 1]), max(p4c_obj$dgram.norm$dgram[scope_idx, 
            1]), mean(p4c_obj$dgram.norm$dgram[scope_idx, 2]), mean(ref_p4c_obj$dgram[scope_idx, 
            2]), "NA")
    }
    
    colnames(slice_df) <- c("start", "end", p4c_obj$track_nm, ref_p4c_obj$track_nm, 
        "mean_scale")
    slice_df$log2_foldChange <- log2(slice_df[, 3]/slice_df[, 4])
    
    # transform to fend space the scales
    slice_df$mean_scale <- round(slice_df$mean_scale * 2)
    
    # calc chisq p-value for region
    scope_idx <- which(temp1$dgram[, 1] >= start & temp1$dgram[, 1] <= end)
    
    N1 <- sum(temp1$dgram[, 2])
    N2 <- sum(temp2$dgram[, 2])
    n1 <- sum(temp1$dgram[scope_idx, 2])
    n2 <- sum(temp2$dgram[scope_idx, 2])
    
    mat <- matrix(c(n1, N1 - n1, n2, N2 - n2), nrow = 2, byrow = F, dimnames = list(c("win", 
        "tot"), c("1", "2")))
    pv <- (chisq.test(mat))$p.value
    
    
    slice_df$p.value <- pv
    return(slice_df)
}
#' Export bedGraph
#' 
#' \code{p4cExportBedGraph} exports the smoothed trend with an adaptative or a fixed smoothing window
#' to a file which can be imported as a 'custom track' to the UCSC genome browser.
#' 
#' @param p4c_obj \code{p4cProfile} object
#' @param filename Filename of the bedGraph.
#' @param trend_scale This parameter controls the smoothing method for the trendline. Either 
#'  'adaptive' (default), or, for fixed window trend, an integer with the desired size 
#'  of the window (in restriction fragments units), put 0 for the raw molecule number.
#' @param min_win_cov If the default smoothing method is used ('adaptive'), this parameter
#'  controls the window size by requiring that no less than \code{min_win_cov} molecules will be included
#'  in that genomic window.
#' @param color Color of the bedGraph plot.
#' 
#' @export
p4cExportBedGraph <- function(p4c_obj, filename, min_win_cov = 50, color = "black")
{
    # Check that the parameters are valid
    if (trend_scale != "adaptive" & !is.numeric(trend_scale))
    {
        stop("illegal trend scale ", trend_scale)
    }
    # Check that the requested trend_scale is part of the dgram_scales
    if (is.numeric(trend_scale) & !trend_scale %in% c(0,p4c_obj$dgram_params$dgram_scales))
    {
        stop("trend_scale was not generated for ", trend_scale, ": set dgram_scales in p4c.conf")
    }
    
    cur_name <- p4c_obj$track_nm
    cur_col <- paste((col2rgb(color)), collapse = ",")

    # Get the data:
    if (is.numeric(trend_scale))
    {
        # access dgram table
        dgram <- p4c_obj$dgram
        # Get the coordinates
        coords <- dgram[, 1]
        bait_x <- p4c_obj$bait$start
        # Define the index of the desired trend_scale in the dgram matrix
        # If needed, update the coordinates
        # As well as a track_name and description
        if (trend_scale == 0)
        {
            trend_scale_idx <- 2
            cur_desc <- sprintf("%s _ UMI-4C track raw values in scope", cur_name)
            track_name <- paste0(cur_name,"_raw")
            message("bedgraph with raw values.")
        } else
        {
            trend_scale_idx <- which(p4c_obj$dgram_params$dgram_scales == trend_scale) + 2
            # Update the coordinates to use geomean coordinates
            coords <- p4cWinGeoMeanCoordinate(coords, trend_scale * 2, bait_x)
            n <- length(dgram[, 1])
            # Restrict to inside values to remove NAs
            coords <- coords[trend_scale:(n - trend_scale)]
            dgram <- dgram[trend_scale:(n - trend_scale), ]
            cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (fixed trend scale: %i fragments)",
                cur_name, trend_scale)
            track_name <- paste0(cur_name,"_smoothed",trend_scale,"frags")
            message("fixed trend scale: ", trend_scale)
        }
        # Extract the trend values
        trend <- dgram[, trend_scale_idx]
    } else if (trend_scale == "adaptive")
    {
        if ("smoothedTrend" %in% names(p4c_obj))
        {
            message("trend is adaptive")
        } else
        {
            message("Will create smoothed trend for track with min_win_cov = ", min_win_cov)
            p4c_obj <- .p4cSmoothedTrend.p4cProfile(p4c_obj, min_win_cov)
        }
        coords <- p4c_obj$smoothedTrend$start
        trend <- p4c_obj$smoothedTrend$trend
        cur_desc <- sprintf("%s _ UMI-4C track smoothed values in scope (min coverage - %i molecules)", 
            cur_name, min_win_cov)
        track_name <- paste0(cur_name,"_smoothedMin",min_win_cov,"mols")
    }
    header_line1 <- sprintf("browser position %s:%i-%i", paste0("chr", gsub("chr", 
        "", p4c_obj$bait$chrom)), round(coords[1]), round(tail(coords, 1)))
    header_line2 <- sprintf("track type=bedGraph name=\"%s\" description=\"%s\" color=%s", 
        track_name, cur_desc, cur_col)
    write(header_line1, filename)
    write(header_line2, filename, append = T)
    
    bed_df <- data.frame(chrom = paste0("chr", gsub("chr", "", p4c_obj$bait$chrom)), 
        start = round(coords), end = round(coords) + 1, 
        value = trend)
    
    write.table(bed_df, filename, append = T, quote = F, col.names = F, row.names = F, sep="\t")
    message("Wrote bedGraph file to ", filename)
}


########################################################################
# Internal functions for p4c package
#######################################################################
.p4cGenerateBaitDgram <- function(p4c_obj) UseMethod(".p4cGenerateBaitDgram")
.p4cGenerateBaitDgram.p4cProfile = function(p4c_obj)
{
    stat_type <- p4c_obj$dgram_params$stat_type
    message("Will use stat_type ", stat_type)
    kill_self_horiz <- getOption("TG3C.kill_self_horiz", 1500)
    
    re_seq <- p4c_obj$dgram_params$re_seq
    map_thresh <- p4c_obj$dgram_params$map_thresh
    min_flen <- p4c_obj$dgram_params$min_flen
    scales <- p4c_obj$dgram_params$dgram_scales
    
    bait_lookup_expansion <- p4c_obj$bait$bait_lookup_expansion
    
    scope_5 <- p4c_obj$scope["scope_5"]
    scope_3 <- p4c_obj$scope["scope_3"]
    
    bait_chrom <- p4c_obj$bait$chrom
    bait_start <- p4c_obj$bait$start
    
    re_map_track = sprintf("redb.fe%s_map", re_seq)
    re_fl_track = sprintf("redb.fe%s_flen", re_seq)
    
    cis_scope = gintervals(bait_chrom, bait_start - scope_5, bait_start + scope_3)
    
    re_scope = gextract(re_map_track, re_fl_track, intervals = cis_scope, iterator = re_map_track, 
        colnames = c("map", "fl"))
    
    re_map = re_scope[re_scope$map > map_thresh & re_scope$fl > min_flen, ]
    
    n = dim(re_map)[1]
    
    re_interv = gintervals.2d(rep(bait_chrom, n), rep(as.numeric(bait_start) - bait_lookup_expansion, 
        n), rep(as.numeric(bait_start) + bait_lookup_expansion, n), re_map$chrom, re_map$start, 
        re_map$end)
    
    # add support for more than one track
    track_nms <- unlist(strsplit(p4c_obj$track_nm, " "))
    vtrack_nms <- paste0("cisprof", 1:length(track_nms))
    for (i in 1:length(track_nms))
    {
        gvtrack.create(vtrack_nms[i], track_nms[i], "weighted.sum")
    }
    expr <- paste(sprintf("ifelse(!is.nan(%s),%s,0)", vtrack_nms, vtrack_nms), collapse = " + ")
    
    prof = gextract(expr, intervals = re_interv, iterator = re_interv, colnames = c("mol"))
    
    prof = prof[abs(prof$start2 - as.numeric(bait_start)) > kill_self_horiz, ]
    bait_idx <- which(abs(prof$start2 - as.numeric(bait_start)) == min(abs(prof$start2 - 
        as.numeric(bait_start))))[1]
    message("Done extract, contacts in window: ", sum(prof$mol))
    
    dgram = matrix(prof$start2)
    n = nrow(prof)
    
    if (stat_type == "linear")
    {
        cprof = cumsum(prof$mol)
        cprof_to_bait <- cumsum(prof$mol[1:bait_idx])
        cprof_from_bait <- cumsum(prof$mol[(bait_idx + 1):n])
        n2 <- length(cprof_from_bait)
    } else if (stat_type == "log")
    {
        cprof_to_bait <- cumsum(log2(1 + prof$mol[1:bait_idx]))
        cprof_from_bait <- cumsum(log2(1 + prof$mol[(bait_idx + 1):n]))
        n2 <- length(cprof_from_bait)
    }
    
    dgram = cbind(dgram, prof$mol)
    
    # check if scope is large enough to run the scales
    if (bait_idx - 2 * max(p4c_obj$dgram_params$dgram_scales) < 2)
    {
        stop("Scope is too small for building the domainogram.\nConsider increasing the genomic scope or reducing the number of scales in p4c_conf.r dgram_scales")
    }
    
    for (d in scales)
    {
        if (stat_type == "linear")
        {
            dgram1 <- c(rep(NA, d), (cprof_to_bait[(2 * d + 1):bait_idx] - cprof_to_bait[1:(bait_idx - 
                2 * d)])/(2 * d), rep(NA, d))
            dgram2 <- c(rep(NA, d), (cprof_from_bait[(2 * d + 1):n2] - cprof_from_bait[1:(n2 - 
                2 * d)])/(2 * d), rep(NA, d))
            dgram <- cbind(dgram, c(dgram1, dgram2))
        } else if (stat_type == "log")
        {
            dgram1 <- c(rep(NA, d), (cprof_to_bait[(2 * d + 1):bait_idx] - cprof_to_bait[1:(bait_idx - 
                2 * d)])/(2 * d), rep(NA, d))
            dgram2 <- c(rep(NA, d), (cprof_from_bait[(2 * d + 1):n2] - cprof_from_bait[1:(n2 - 
                2 * d)])/(2 * d), rep(NA, d))
            dgram <- cbind(dgram, c(dgram1, dgram2))
        }
    }
    
    n = dim(dgram)[1]
    p4c_obj$dgram <- (dgram[10:(n - 10), ])
    p4c_obj
}


#### Internal function for returning the bin mols sum for normalization purpose

.p4cGenCovNormVector <- function(p4c_obj, norm_bins = 10^(3:6), r_expand = 1.2, post_smooth_win = 50, 
    stat_type = "log")
    {
    bins <- sort(norm_bins)
    
    p4c_bait_start <- p4c_obj$bait$start
    chrom <- p4c_obj$bait$chrom
    scope <- p4c_obj$scope
    stat_type <- p4c_obj$dgram_params$stat_type
    track_nms <- unlist(strsplit(p4c_obj$track_nm, " "))
    
    # Check that the scope is shorter than max norm_bin
    if (max(scope) > max(bins))
    {
        message(sprintf("Error: Max normalizing bin (%i) is smaller than max scope (%i)", 
            max(bins), max(scope)))
        message("Increase maximal value of TG3C.norm_bins in p4c.conf")
        stop()
    }
    
    expanded_obj <- p4cNewProfile(track_nms, chrom, p4c_bait_start, max(bins), max(bins), 
        stat_type = stat_type)
    
    dgram <- p4c_obj$dgram
    
    expanded_dgram <- expanded_obj$dgram
    mols = expanded_dgram[, 2]
    
    coords = expanded_dgram[, 1]
    dlt = abs(coords - p4c_bait_start)
    
    sum_vec <- rep(NA, length(coords))
    
    f_bin = dlt < bins[1]
    tot_in_bin <- sum(log2(1 + mols[f_bin]))
    if (stat_type == "linear") 
        tot_in_bin <- sum(mols[f_bin])
    sum_vec[f_bin] = tot_in_bin
    
    for (i in 2:(length(bins)))
    {
        f_bin_exp = dlt < bins[i] * r_expand & dlt >= bins[i - 1]/r_expand
        f_bin = dlt < bins[i] & dlt >= bins[i - 1]
        if (stat_type == "log"){
            tot_in_bin = sum(log2(1 + mols[f_bin_exp]))
        } else if (stat_type == "linear"){ 
            tot_in_bin <- sum(mols[f_bin_exp])
        }
        sum_vec[f_bin] = tot_in_bin
        message("bin ", bins[i], " tot ", tot_in_bin)
    }

    f_bin = dlt >= bins[length(bins)]
    if (stat_type == "log"){
        tot_in_bin <- sum(log2(1 + mols[f_bin]))
    } else if (stat_type == "linear"){
        tot_in_bin <- sum(mols[f_bin])
    }
    sum_vec[f_bin] = tot_in_bin
    
    sum_vec = zoo::rollmean(sum_vec, post_smooth_win, fill = c(sum_vec[1], NA, sum_vec[length(sum_vec)]))
    # return only scope values
    f_coords_in_scope <- coords %in% dgram[, 1]
    
    
    return(sum_vec[f_coords_in_scope])
}



.p4cNormDgram <- function(p4c_obj, ref_p4c_obj, norm_bins) UseMethod(".p4cNormDgram")

.p4cNormDgram.p4cProfile = function(p4c_obj, ref_p4c_obj, norm_bins = getOption("TG3C.norm_bins"))
{
    if (p4c_obj$scope["scope_5"] != ref_p4c_obj$scope["scope_5"] | p4c_obj$scope["scope_3"] != 
        ref_p4c_obj$scope["scope_3"])
        {
        stop("scope of ref p4c object is different")
    }
    
    if (p4c_obj$bait$chrom != ref_p4c_obj$bait$chrom | p4c_obj$bait$start != ref_p4c_obj$bait$start)
    {
        stop("bait of ref p4c object is different")
    }
    
    if (is.null(p4c_obj[["dgram"]]))
    {
        stop("no dgram matrix for p4c_obj")
    }
    
    if (is.null(ref_p4c_obj[["dgram"]]))
    {
        stop("no dgram matrix for ref_p4c_obj")
    }
    
    bait_start <- p4c_obj$bait$start
    stat_type <- p4c_obj$dgram_params$stat_type
    dgram1 <- p4c_obj$dgram
    dgram2 <- ref_p4c_obj$dgram
    
    tot1 <- .p4cGenCovNormVector(p4c_obj, stat_type = stat_type, norm_bins = norm_bins)
    tot2 <- .p4cGenCovNormVector(ref_p4c_obj, stat_type = stat_type, norm_bins = norm_bins)
    
    mols1 <- dgram1[, 2]
    mols2 <- dgram2[, 2]
    
    norm1 <- tot2/tot1
    
    dgram1[, -1] <- dgram1[, -1] * norm1
    
    p4c_obj$dgram.norm <- list(dgram = dgram1, factors = norm1)
    p4c_obj$normalized_by <- ref_p4c_obj$track_nm
    return(p4c_obj)
}

# Adaptive smoothing function for a single profile
.p4cSmoothedTrend <- function(p4c_obj, min_win_cov, min_res) UseMethod(".p4cSmoothedTrend")
.p4cSmoothedTrend.p4cProfile <- function(p4c_obj, min_win_cov = 30, min_res = 1)
{
    dgram <- p4c_obj$dgram
    # passing over all values in the matrix and looking for the min values which pass
    # threshold filling the smoothing vector with this values.
    
    dgram[is.na(dgram)] <- 0
    bait_x <- p4c_obj$bait$start
    vec1 <- rep(NA, nrow(dgram))
    scale <- rep(NA, nrow(dgram))
    coord <- rep(NA, nrow(dgram))
    base_coord <- dgram[, 1]
    mean_coord <- numeric()
    base_scales <- p4c_obj$dgram_params$dgram_scales
    
    for (i in (2 + min_res):ncol(dgram))
    {
        cur_scale <- base_scales[i - 2]
        f <- is.na(vec1) & (dgram[, i] * cur_scale * 2 > min_win_cov)
        vec1[f] <- dgram[f, i]
        mean_coords <- .p4cWinGeoMeanCoordinate(base_coord, cur_scale * 2, bait_x)
        coord[f] <- mean_coords[f]
        scale[f] <- cur_scale
    }
    
    vec1[is.na(vec1)] <- NA
    coord[is.na(coord)] <- base_coord[is.na(coord)]
    p4c_obj$smoothedTrend <- list(start = sort(coord), trend = vec1, scale = scale)
    return(p4c_obj)
}

# Adaptive smoothing function for a comparative profile
.p4cSmoothedTrendComp <- function(p4c_obj1, ref_p4c_obj, min_win_cov, min_res) UseMethod(".p4cSmoothedTrendComp")
.p4cSmoothedTrendComp.p4cProfile <- function(p4c_obj1, ref_p4c_obj, min_win_cov = 30, 
    min_res = 1)
    {
    # First do adapted smoothing, and after normalize by read count factor
    # (norm_factors)
    dgram1.norm <- p4c_obj1$dgram.norm$dgram
    norm_factors <- p4c_obj1$dgram.norm$factors
    dgram1 <- p4c_obj1$dgram
    dgram2 <- ref_p4c_obj$dgram
    
    dgram1[is.na(dgram1)] <- 0
    dgram1.norm[is.na(dgram1.norm)] <- 0
    dgram2[is.na(dgram2)] <- 0
    
    if (is.null(p4c_obj1[["dgram.norm"]]))
    {
        stop("no normalized matrix, use .p4cNormDgram")
    }
    
    if (p4c_obj1$normalized_by != ref_p4c_obj$track_nm)
    {
        stop("p4c_obj1 is not normalized to ref_p4c_obj")
    }
    
    ## create adaptive trend
    bait_x <- p4c_obj1$bait$start
    vec1 <- rep(NA, nrow(dgram1))
    vec2 <- rep(NA, nrow(dgram2))
    coord <- rep(NA, nrow(dgram1))
    base_coord <- dgram1[, 1]
    scale <- rep(NA, nrow(dgram1))
    mean_coord <- numeric()
    base_scales <- p4c_obj1$dgram_params$dgram_scales
    
    # passing over all values in the matrix and looking for the min values which pass
    # threshold filling the smoothing vector with this values.
    message("Will attempt getting coverage for ", ncol(dgram1), " scales")
    for (i in (2 + min_res):ncol(dgram1))
    {
        cur_scale <- base_scales[i - 2]
        f <- is.na(vec1) & 
            (dgram1[, i] * cur_scale * 2 > min_win_cov & 
            dgram2[, i] * cur_scale * 2 > min_win_cov)
        
        vec1[f] <- dgram1[f, i]
        vec2[f] <- dgram2[f, i]
        mean_coords <- .p4cWinGeoMeanCoordinate(base_coord, cur_scale * 2, bait_x)
        coord[f] <- mean_coords[f]
        
        scale[f] <- cur_scale
    }
    
    vec1.norm <- vec1 * norm_factors  #normalize trend of obj1
    coord[is.na(coord)] <- base_coord[is.na(coord)]
    trend_mat <- data.frame(start = sort(coord), p4c_obj1 = vec1.norm, ref_p4c_obj = vec2, 
        scale = scale)
    dlt_dgram <- log2(1 + dgram1.norm[, -1]) - log2(1 + dgram2[, -1])
    
    return(list(trend_mat = trend_mat, dlt_dgram = dlt_dgram))
}



# Calc geometric mean of coordinates relative to distance from bait
.p4cWinGeoMeanCoordinate <- function(coords, scale, bait_x)
{
    bait_idx <- which(abs(coords - as.numeric(bait_x)) == min(abs(coords - as.numeric(bait_x))))[1]
    offsets <- abs(coords - coords[bait_idx])
    offsets[bait_idx] <- 1  #avoiding plugging 0s in the log
    
    
    mean_offsets1 <- exp(zoo::rollmean(log(offsets[1:bait_idx]), scale, fill = NA))
    mean_offsets2 <- exp(zoo::rollmean(log(offsets[(bait_idx + 1):length(offsets)]), scale, 
        fill = NA))
    
    # deal with NAs near the bait
    near_bait1 <- sapply((bait_idx - scale/2 + 1):bait_idx, function(i) exp(mean(log(offsets[i:bait_idx]))))
    near_bait2 <- sapply((bait_idx + 1):(bait_idx + scale/2 - 1), function(i) exp(mean(log(offsets[(bait_idx + 
        1):i]))))
    
    mean_offsets1[(bait_idx - scale/2 + 1):bait_idx] <- near_bait1
    mean_offsets2[1:(scale/2 - 1)] <- near_bait2
    
    mean_offsets1 <- mean_offsets1 * -1
    mean_offsets1[bait_idx] <- 0
    
    mean_offsets <- c(mean_offsets1, mean_offsets2)
    m_coords <- coords[bait_idx] + mean_offsets
    
    return(m_coords)
}

.p4cSetCurBait <- function(p4c_obj, bait_chrom, bait_start) UseMethod(".p4cSetCurBait")
.p4cSetCurBait.p4cProfile <- function(p4c_obj, bait_chrom = NULL, bait_start = NULL, 
    bait_lookup_expansion = getOption("TG3C.bait_lookup_expansion"))
    {
    if (is.null(bait_chrom) | is.null(bait_start))
    {
        message("no new bait info, current bait:")
        return(p4c_obj$bait)
    }
    
    p4c_obj$bait <- list(chrom = bait_chrom, start = as.numeric(bait_start), bait_lookup_expansion = as.numeric(bait_lookup_expansion))
    if (!is.null(p4c_obj[["dgram"]]))
    {
        message("updating Dgram matrix")
        p4c_obj <- .p4cGenerateBaitDgram(p4c_obj)
    }
    
    if (!is.null(p4c_obj[["dgram.norm"]]))
    {
        warning("changed current bait, the normalized dgram is no longer valid!\n
            use .p4cNormDgram again")
    }
    return(p4c_obj)
}

.p4cSetCurScope <- function(p4c_obj, scope_5, scope_3) UseMethod(".p4cSetCurScope")
.p4cSetCurScope.p4cProfile <- function(p4c_obj, scope_5 = NULL, scope_3 = NULL)
{
    if (missing(scope_5) | missing(scope_3))
    {
        message("no scope supplied, current scope:")
        return(p4c_obj$scope)
    }
    
    p4c_obj$scope <- c(scope_5 = scope_5, scope_3 = scope_3)
    if (!is.null(p4c_obj[["dgram"]]))
    {
        message("updating Dgram matrix")
        p4c_obj <- .p4cGenerateBaitDgram(p4c_obj)
    }
    
    if (!is.null(p4c_obj[["dgram.norm"]]))
    {
        warning("changed current bait, the normalized dgram is no longer valid!\n
            use .p4cNormDgram again")
    }
    return(p4c_obj)
}

.p4cSetDgramParams <- function(p4c_obj, re_seq, map_thresh, min_flen, scales, stat_type) UseMethod(".p4cSetDgramParams")
.p4cSetDgramParams.p4cProfile <- function(p4c_obj, re_seq = getOption("TG3C.RE_seq"), 
    map_thresh = getOption("TG3C.map_thresh"), min_flen = getOption("TG3C.min_flen"), 
    scales = getOption("TG3C.dgram_scales"), stat_type = getOption("TG3C.stat_type"))
{
    l <- list(re_seq = re_seq, map_thresh = as.numeric(map_thresh), min_flen = as.numeric(min_flen), 
        dgram_scales = scales, stat_type = stat_type)
    p4c_obj$dgram_params <- l
    return(p4c_obj)
}

.p4cSetPngParams <- function(p4c_obj, png_w, png_h, png_dir, png_res) UseMethod(".p4cSetPngParams")
.p4cSetPngParams.p4cProfile <- function(p4c_obj, png_w = getOption("TG3C.png_w"), png_h = getOption("TG3C.png_h"), 
    png_dir = getOption("TG3C.png_dir"), png_res = getOption("TG3C.png_res"))
{
    l <- list(png_w = as.numeric(png_w), png_h = as.numeric(png_h), png_dir = png_dir, 
        png_res = as.numeric(png_res))
    p4c_obj$figure_params <- l
    return(p4c_obj)
}

.p4cSetDgramGraphics <- function(p4c_obj, zlim, shades) UseMethod(".p4cSetDgramGraphics")
.p4cSetDgramGraphics.p4cProfile <- function(p4c_obj, zlim = getOption("TG3C.p4c_dgram_zlim"), 
    shades = getOption("TG3C.p4c_dgram_shades"))
{
    p4c_obj$graphics_params <- list(zlim = zlim, shades = shades)
    return(p4c_obj)
}

# color scale bar function:
.image.scale <- function(z, zlim, col = heat.colors(12), breaks, horiz = TRUE, ylim = NULL, 
    xlim = NULL, yaxt = "n", ...)
{
    if (!missing(breaks)) {
        if (length(breaks) != (length(col) + 1))
        {
            stop("must have one more break than colour")
        }
    }
    
    if (missing(breaks) & !missing(zlim)) {
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
    }
    
    if (missing(breaks) & missing(zlim)) {
        zlim <- range(z, na.rm = TRUE)
        zlim[2] <- zlim[2] + c(zlim[2] - zlim[1]) * (0.001)
        zlim[1] <- zlim[1] - c(zlim[2] - zlim[1]) * (0.001)
        breaks <- seq(zlim[1], zlim[2], length.out = (length(col) + 1))
    }
    
    poly <- vector(mode = "list", length(col))
    
    for (i in seq(poly)) {
        poly[[i]] <- c(breaks[i], breaks[i + 1], breaks[i + 1], breaks[i])
    }
    
    xaxt <- ifelse(horiz, "s", "n")
    yaxt <- ifelse(horiz, "n", "s")
    
    if (horiz) {
        YLIM <- c(0, 1)
        XLIM <- range(breaks)
    }
    
    if (!horiz) {
        YLIM <- range(breaks)
        XLIM <- c(0, 1)
    }
    
    if (missing(xlim)){ 
        xlim = XLIM
    }
    if (missing(ylim)){
        ylim = YLIM
    }
    
    plot(1, 1, t = "n", ylim = ylim, xlim = xlim, xaxt = xaxt, yaxt = yaxt, xaxs = "i", 
        yaxs = "i", ...)
    for (i in seq(poly)) {
        if (horiz) {
            polygon(poly[[i]], c(0, 0, 1, 1), col = col[i], border = NA)
        }
        if (!horiz) {
            polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
        }
    }
}

