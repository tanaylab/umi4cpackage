# A group of routines for 3C/4C tracks

#' Build the required genomic datasets.
#' 
#' Downloads the required datasets for the package: trackdb and redb. These datasets
#' include a genomic database in a binary format and restriction sites tables in text format.
#' 
#' @param output_dir Output directory for the dataset.
#' @param genome For which genome build should the data be download. Current options available are:
#' 'hg19', 'mm9', 'mm10'.
#' @param reseq Restriction enzye database. Currently only DpnII (\strong{GATC}) is available.
#' 
#' @export 
p4cBuildRequirements <- function(output_dir, genome="hg19", reseq="GATC", 
                                 dataset_keys="http://www.wisdom.weizmann.ac.il/~omersch/dataset_keys.txt")
{
    datasets <- read.table(url(dataset_keys), header=TRUE, stringsAsFactors=FALSE)
    tf <- tempfile()
    dir.create(output_dir, recursive=TRUE)
    
    loc <- datasets[ datasets$genome == genome & datasets$reseq == reseq, 'loc']
    if (length(loc) != 1) 
    {
        stop(sprintf("Couldn't find datasets for genome %s and reseq %s", genome, reseq))
    }
    
    message("Will begin downloading datasets to ", tf)
    ret <- download.file(loc, tf, mode='wb')
    if (ret != 0)
    {
        stop("Couldn't download file from ", loc)
    }
    
    untar(tf, exdir=output_dir, compressed=TRUE, verbose=TRUE)
    message("Done writing genomic dataset files to ", output_dir)
}

#' Dump config files templates.
#' 
#' Dump templates of config files required by the package.
#' 
#' @param conf_dir Directory to dump to files to.
#' 
#' @export 
p4cDumpConfFiles <- function(conf_dir)
{
    conf_files <- dir(system.file("conf", package='umi4cPackage'), full.names=T)
    dir.create(conf_dir, recursive=T, showWarnings=FALSE)
    ret <- file.copy(conf_files, conf_dir, recursive=T, overwrite=FALSE)
    if (!all(ret)) 
    {
        warning("couldn't dump conf files to ", conf_dir, "\n  Perhaps they're already there? ")    
    } else
    {
        message("Dumped conf (config) files to: ", conf_dir)
    }
    message("\nMake sure you fill all the conf files before continuing!\n")  
} 

#' Load all the config files required for the package.
#' 
#' Loads the config files which were filled by the user. Accepts either directory path
#' from \code{conf_dir}, or individual files supplied by the other parameters.
#' 
#' @param conf_dir A directory with all the defined conf files. Should have the same
#' files as those that were exported with \code{p4cDumpConfFiles}.
#' @param paths_conf Parameter file like \emph{paths.conf}
#' @param pipeline_conf Parameter file like \emph{pipeline.conf}
#' @param p4c_conf Parameter file like \emph{p4c.conf}
#' 
#' @export 
p4cLoadConfFiles <- function(conf_dir, 
                             paths_conf=file.path(conf_dir, "paths.conf"),
                             pipeline_conf=file.path(conf_dir, "pipeline.conf"),
                             p4c_conf=file.path(conf_dir, "p4c.conf"))
{
    op <- options()
    all_paths <- list(
        "TG3C.pipeline_conf_path"=pipeline_conf,
        "TG3C.paths_conf_path"=paths_conf,
        "TG3C.p4c_conf_path"=p4c_conf)
    
    op.p4c <- init_params(paths_conf)
    op.p4c <- c(op.p4c, init_params(pipeline_conf), init_params(p4c_conf))
    op.p4c <- c(op.p4c, all_paths)
    
    # Check if some params where saved alread
    setted <- names(op) %in% names(op.p4c)
    if (any(setted))
    {
        
        message(sprintf("%i parameters were already defined in this sessions and they will overwritten!", 
                        sum(setted)))
    }
    
    options(op.p4c)
    
    .checkConf()
    
    gsetroot(getOption("TG3C.trackdb"), rescan = TRUE)
    
    invisible()
}


# Check conf files were loaded.
.checkConf <- function(additional_params = NULL)
{
    params <- c("TG3C.trackdb", "TG3C.kill_self_horiz", "TG3C.switch_ratio")
    if (!is.null(additional_params)) 
    {
        params <- c(params, additional_params)
    }
    
    if (!all( params %in% names(options()) )) 
    {
        stop("Missing conf parameters! Have you loaded the files with 4cLoadConfFiles?")
    }
}

#' Import fastq files to a genomic track.
#' 
#' Generate a genomic track from fastqs listed in the samples.conf file.
#' 
#' @param sample_ids Numerical samples IDs from the baits.txt table.
#' @param track_desc Description attribute. Optional.
#' @param overwrite.if.exists Force overwriting on older tracks.
#' @param verbose Logical. Controls verbosity. 
#' @param groot Optional. Controls the trackdb to which tracks are imported.
#' 
#' @examples
#' \donttest{
#' # Load config files
#' p4cLoadConfFiles("conf/")
#' 
#' # Import 3 tracks
#' p4cCreate4CseqTrack(1:3)
#' }
#' 
#' @export 
p4cCreate4CseqTrack = function(sample_ids = NULL, track_desc = "4C track",  verbose = TRUE,  
                               overwrite.if.exist = FALSE, groot = NULL)
{
    .checkConf()
    options(gparam.type = "string")
    Sys.setenv(PERL_BADLANG = 0)
    pipeline_conf_path <- getOption("TG3C.pipeline_conf_path")
    paths_conf_path    <- getOption("TG3C.paths_conf_path")
    perl_scripts_path <- system.file("perl", package="umi4cPackage")
    
    
    # read table data from tables
    baits_tab_nm <- getOption("TG3C.baits_tab")
    samples_tab_nm <- getOption("TG3C.samples_tab")
    
    baits_tab <- read.table(baits_tab_nm, header = TRUE, stringsAsFactors=FALSE)
    if (anyDuplicated(baits_tab$Bait_ID)) 
    {
        stop("Duplicated values in baits table - make sure not duplicated values in bait_ID")
    }
    samples_tab <- read.table(samples_tab_nm, header = TRUE, stringsAsFactors=FALSE)
    if (anyDuplicated(samples_tab$Sample_ID)) 
    {
        stop("Duplicated values in samples table - make sure not duplicated values in sample_ID")
    }
    
    
    trackdb_path = getOption("TG3C.trackdb")
    imp_3C_pipe_pl = getOption("TG3C.imp_3C_pipe_pl")
    redb = getOption("TG3C.redb")
    re_seq = getOption("TG3C.RE_seq")
    fends_fn = sprintf("%s/%s.fends", redb, re_seq)
    base_track_nm = getOption("TG3C.base_track_name")
    
    if (is.null(groot))
    {
        gsetroot(trackdb_path, rescan = TRUE)
        message("Will write tracks to ", trackdb_path)
    }
    
    work_dir <- getOption("TG3C.workdir")
    
    message("Work dir is ", work_dir)
    message("\n\nImport 4C from fastq\n")
    
    
    if (!file.exists(fends_fn))
    {
        message("ERROR: fends table is missing from ", fends_fn)
    }
    
    adj_list = list()
    if (is.null(sample_ids)){
        sample_ids <- samples_tab$Sample_ID
    }
    for (foc_ndx in sample_ids)
    {
        row_idx <- which(samples_tab$Sample_ID == foc_ndx)
        
        # Test existence
        if (length(row_idx) == 0) 
        {
            stop("No Sample_ID ", foc_ndx, ", check samples table!")
        }
        if (length(row_idx) > 1)
        {
            stop("Sample_ID in samples table are not unique! (In Sample_ID ", foc_ndx, " )")
        }
        
        message("working on sample ID ", foc_ndx)
        exp_nm <- samples_tab$Experiment_name[row_idx]
        sample_nm <- samples_tab$Sample_name[row_idx]
        baits_idx <- as.numeric(strsplit(as.character(samples_tab$Bait_IDs[row_idx]), 
                                         split = ",")[[1]])
        
        foc_ndx_track_name = sprintf("%s_%s_%s", base_track_nm, exp_nm, sample_nm)
        
        # Check existence  
        existing_tracks <- 0
        for (bait_id in baits_idx)
        {
            bait_row_idx <- which(baits_tab$Bait_ID == bait_id)
            if (length(bait_row_idx) == 0) 
            {
                stop("No Bait_ID ", bait_id, " for Sample_ID ", foc_ndx, "\n  Check baits/samples tables!")
            }
            if (length(bait_row_idx) > 1)
            {
                stop("Bait_ID in baits table are not unique! (In Bait_ID ", bait_id, " )")
            }
            
            bait_name <- baits_tab$Bait_name[bait_row_idx]
            foc_ndx_track_name_bait <- paste0(foc_ndx_track_name, "_", bait_name)
            if (overwrite.if.exist & gtrack.exists(foc_ndx_track_name_bait))
            {
                message("track ", foc_ndx_track_name_bait, " already exist, overwriting it...\n")
                gtrack.rm(foc_ndx_track_name_bait, T)
            }
            if (gtrack.exists(foc_ndx_track_name_bait) & !overwrite.if.exist)
            {
                message("Skipping processing of track ", foc_ndx_track_name_bait, ": already exists.\n")
                existing_tracks <- existing_tracks + 1
                next
            } 
        }
        if (existing_tracks == length(baits_idx)) 
        {
            message('All tracks from sample ID ', foc_ndx, ' already imported, skipping sample ID...\n')
            next
        }
        
        cmd = sprintf("perl %s/%s @%s @%s -TG3C.sample_id %s -TG3C.perl_scripts_path %s\n", 
                      perl_scripts_path, 
                      imp_3C_pipe_pl, 
                      pipeline_conf_path, 
                      paths_conf_path, 
                      foc_ndx, 
                      perl_scripts_path)
        if (verbose)
        {
            message(cmd, "\n")
        }
        ret <- system(cmd)
        if (ret)
        {
            stop("An error occured while running the pipeline on index ", foc_ndx, 
                 " ")
        }
        
        status_fn = sprintf("%s/%s.%s/done_ok", work_dir, exp_nm, sample_nm)
        if (!file.exists(status_fn))
        {
            message("ERROR: failed in index ", foc_ndx)
            message(sprintf("See more logs in workdir: %s/%s.%s/logs ", work_dir, exp_nm, sample_nm))
            return(NA)
        }
        
        
        # Genereate a track for each bait
        for (bait_id in baits_idx)
        {
            bait_row_idx <- which(baits_tab$Bait_ID == bait_id)
            bait_name <- baits_tab$Bait_name[bait_row_idx]
            contacts_fn <- sprintf("%s/%s.%s/adj.%s", work_dir, exp_nm, sample_nm, 
                                   bait_name)
            foc_ndx_track_name_bait <- paste0(foc_ndx_track_name, "_", bait_name)
            if (file.exists(contacts_fn)) {
                message("Will import ", foc_ndx_track_name_bait, " from ", contacts_fn, 
                        " fends at ", fends_fn)
                gtrack.2d.import_contacts(foc_ndx_track_name_bait, description = track_desc, 
                                          contacts = contacts_fn, fends = fends_fn)
                bait_chr <- baits_tab$Bait_chr[bait_row_idx]
                bait_coord <- baits_tab$Bait_coord[bait_row_idx]
                gtrack.attr.set(foc_ndx_track_name_bait, attr = "Bait_name", value = bait_name)
                gtrack.attr.set(foc_ndx_track_name_bait, attr = "Bait_chr", value = bait_chr)
                gtrack.attr.set(foc_ndx_track_name_bait, attr = "Bait_coord", value = bait_coord)
                message(foc_ndx_track_name_bait, " was imported succsesfully! ")
            } else {
                message("Will not import ", foc_ndx_track_name_bait, " as ", contacts_fn, 
                        " does not exists")
            }
        }
    }
}

#' Create a new redb for UMI-4C
#' 
#' Generate a new database for a restriction enzyme in the redb.conf file.
#' 
#' @param re_seq sequence recognition (e.g. 'GATC' for DpnII)
#' @param redb_params_fn File with the redb parameters ('redb.conf')
#' 
#' @export 
gtrack.create_redb_tracks = function(re_seq, redb_params_fn, verbose = TRUE)
{
    options(gparam.type = "string")
    Sys.setenv(PERL_BADLANG = 0)
    # Set the parameters in the file:
    op.new = init_params(redb_params_fn)
    options(op.new)
    
    # Get the path of the perl scripts
    gen_re_frags_pl <- grep("build_re_db.pl",
                            list.files(system.file("perl",
                                                   package = "umi4cPackage"),
                                       recursive=T, full.names=T),
                            value=T)
    
    re_frags_to_fends_pl <- grep("re_frags_to_fends_tab.pl",
                                 list.files(system.file("perl",
                                                        package = "umi4cPackage"),
                                            recursive=T, full.names=T),
                                 value=T)
    
    mapab_track = getOption("TG3C.mapab_track", NA)
    re_workdir = getOption("TG3C.re_workdir", NA)
    chrom_key = getOption("TG3C.chrom_seq_key", NA)
    
    # Check existance
    all_set <- all(!is.na(c(mapab_track, 
                            re_workdir, chrom_key))) &
        length(gen_re_frags_pl) == 1 &
        length(re_frags_to_fends_pl) == 1
    
    if (!all_set)
    {
        stop("ERROR: missing configuration options\n")
    }
    # Update the database
    gdb.reload()
    # create the frags tables
    if (verbose){
        message("constructing table from genomic sequence\n")
    }
    # Create the directory if it does not exists
    dir.create(re_workdir, recursive = T, showWarnings = F)
    cmd = sprintf("perl %s %s %s %s", gen_re_frags_pl, re_seq, chrom_key, re_workdir)
    system(cmd)
    
    
    # import the frags as sparse tracks
    if (verbose){
        message("importing fragments as sparse tracks")
    }
    
    # create redb dir
    gdir.create("redb", showWarnings = FALSE)
    
    flen_track = sprintf("redb.%s_flen", re_seq)
    gc_track = sprintf("redb.%s_gc", re_seq)
    fragmap_track = sprintf("redb.%s_map", re_seq)
    flen_file = sprintf("%s/%s_flen", re_workdir, re_seq)
    gc_file = sprintf("%s/%s_gc", re_workdir, re_seq)
    fragmap_file = sprintf("%s/%s_map", re_workdir, re_seq)
    if (verbose){
        message("Importing tracks flen_track ", flen_track, " gc_track ", gc_track)
    }
    gtrack.import(flen_track, paste(re_seq, "fragment length"), flen_file, 0)
    gtrack.import(gc_track, paste(re_seq, "gc content"), gc_file, 0)
    # annotate mapability
    if (verbose){
        message("generating mapability data for fraqments, using ", mapab_track)
    }
    gtrack.create(fragmap_track, paste(re_seq, "fragment mapability"),
                  as.character(mapab_track), iterator = flen_track)
    # dump mapability frags to text
    if (verbose){
        message("writing mapability text file to re workdir")
    }
    gextract(fragmap_track, .misha$ALLGENOME, iterator = fragmap_track, file = fragmap_file)
    # separate frags to fends table
    message("writing fend tables from frag tables")
    cmd = sprintf("perl %s %s %s", re_frags_to_fends_pl, re_workdir, re_seq)
    system(cmd)
    
    # import fend tables
    fe_flen_track = sprintf("redb.fe%s_flen", re_seq)
    fe_gc_track = sprintf("redb.fe%s_gc", re_seq)
    fe_fragmap_track = sprintf("redb.fe%s_map", re_seq)
    fe_flen_file = sprintf("%s/fe%s_flen", re_workdir, re_seq)
    fe_gc_file = sprintf("%s/fe%s_gc", re_workdir, re_seq)
    fe_fragmap_file = sprintf("%s/fe%s_map", re_workdir, re_seq)
    # I noticed that when I am using contigs 
    # they may have the re_seq at the beginning: for example in galGal5:
    # >chr17_NT_462588v1_random
    # GGATCCAACAGccaaacagagcaaaaccaaaggaaagtgaggaaaaataa
    # Then 2-3 is both the start and the end
    # of the first fragment. This induce duplicated lines in feFILES.
    nb.duplicates <- system(paste("cat", fe_flen_file,
                                  "| uniq -d | wc -l"),
                            intern = TRUE)
    if (nb.duplicates > 0){
        # If there are duplicates new files will be created ended by
        # wd (without duplicates)
        temp.fe_flen_file <- sprintf("%s/fe%s_flen.wd", re_workdir, re_seq)
        system(paste("cat", fe_flen_file, "| uniq >", temp.fe_flen_file))
        fe_flen_file <- temp.fe_flen_file
        
        temp.fe_gc_file <- sprintf("%s/fe%s_gc.wd", re_workdir, re_seq)
        system(paste("cat",fe_gc_file, "| uniq >", temp.fe_gc_file))
        fe_gc_file <- temp.fe_gc_file
        
        temp.fe_fragmap_file<-sprintf("%s/fe%s_map.wd", re_workdir, re_seq)
        system(paste("cat", fe_fragmap_file, "| uniq >", temp.fe_fragmap_file))
        fe_fragmap_file <- temp.fe_fragmap_file
    }
    if (verbose){
        message(sprintf("Importng fend tracks: %s, %s, %s", fe_flen_track, fe_gc_track, 
                        fe_fragmap_track))
    }
    gtrack.import(fe_flen_track, paste(re_seq, "fragment end length"), fe_flen_file, 
                  0)
    gtrack.import(fe_gc_track, paste(re_seq, "fragment end gc"), fe_gc_file, 0)
    gtrack.import(fe_fragmap_track, paste(re_seq, "fragment end mapability"), fe_fragmap_file, 
                  0)
}

