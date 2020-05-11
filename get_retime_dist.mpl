use CadConfig ;
use ChipVars ;
use Parallel::Loops ;
use Graph::Directed;
use Graph::Traversal::BFS;

# global varables for retiming attributions

our %M_noscan_port_mapping ;
our %M_routeRules ;
our %M_routeRules_pipe ;
our %M_all_neighbor_pars ;
our %M_all_thr_pars ;

our %chiplet_uc ;
our $log_file ;
our $debug ;
our $s_type ;

####define the abut check threshold value, less than 10 treat as not abut partition ######
#set_abutment_dist 10; ## PR gurad band is 10, user could modify it; 

Tsub generate_report => << 'END';
    DESC {
        generate retime detour report, -histogram_min is the min cutoff for distance summary.
        Example:
            MENDER > generate_report NV_gaa_g1.2018Dec18_23_48_fullSM.rep.flat -ultra 1 -skip_clock "jtag*"
            MENDER > generate_report NV_gaa_g1.2018Dec18_23_48_fullSM.rep.flat -histogram_min 700 -histogram_max 2000 -retime_distance 700 -ultra 1
            MENDER > generate_report nv_top.2019Jul30_20.rep.flat -inter_path_only 1
    }
    ARGS {
        -histogram_min:$histogram_min           ## min cutoff for distance historgram, use 700 as default
        -histogram_max:$histogram_max           ## max cutoff for distance historgram, use 2000 as default
        -histogram_step:$histogram_step         ## min cutoff for distance historgram, use 100 as scale step
        -retime_distance:$retime_distance       ## per retime stage distance, use 1200 as default for tsmc16ff process, use 900 as default for 7nm process
        -skip_clock:@skip_clock                 ## skip these clock when dump report, wildchard surported
        -only_clock:@only_clock                 ## only care these clocks when dump report,  wildchard surported
        -use_multi_thread_num:$multi_thread_num ## use multi-thread, default 4 thread 
        -inter_chiplet_only:$opt_inter_only     ## analysis inter-chiplet paths only for nv_top
        -ultra:$opt_ultra                       ## by default ultra mode is off , when -ultra 1, flow enables calculate MCP and mapping hier-unit
        -ipo_dir:$ipo_dir                       ## to load the def/region files @ specified ipo_dir
        -legalize_region_cells                  ## to legalize the cells in regions, skipped the legalization by default 
        -no_clr                                 ## no need to reload the violation files. not default. 
	    -dump_mcp_tcl                           ## to dump the mcp tcl file, default not. 
        -debug                                  ## to print the debug infomations
        -mail:$mail_id                               ## to send a mail to user when done.
        $vio_file 
    }

    ###############################
    ## user defined attributions ##
    ###############################

    define_vio_attr (-class => "path", -is_num => 1, -attr => 'man_distance') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'detour_ratio') ;
    define_vio_attr (-class => "path", -is_num => 0, -attr => 'feed_pars') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'feed_pars_num') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'par_num') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'real_distance') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'ideal_distance') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'is_mcp') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'mcp_setup_num') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'bin_man_distance') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'bin_ideal_distance') ;
    define_vio_attr (-class => "path", -is_num => 1, -attr => 'bin_real_distance') ;
    define_vio_attr (-class => "path", -is_num => 0, -attr => 'sig_name') ;


    ##############################
    ## Options and Environments ##
    ##############################

    my $rep_name = $vio_file ;

    if ($rep_name =~ /\.gz$/) {
        $rep_name =~ s/\.gz$// ;
    }


    $log_file = "$rep_name.gen_retime_rep.log" ; 
    open_log ("$log_file") ;
    
    print "Log File : $log_file\n" ;
    if (defined $opt_debug) {
        $debug = $opt_debug ;
    }
     
    my $top      = get_top;
    my $hist_min = $histogram_min;
    set_chip_top $top;
    $s_type = session_type() ;
    chomp $s_type ;

    lprint "# This is $s_type session. \n\n";

    lprint "Loading def/region files...\n" ;
    if (defined $ipo_dir && defined $opt_legalize_region_cells) {
        load_def_region_files(-ipo_dir => $ipo_dir, -legalize_region) ;
    } elsif (defined $ipo_dir) {
        load_def_region_files(-ipo_dir => $ipo_dir) ;
    } elsif (defined $opt_legalize_region_cells) {
        load_def_region_files(-legalize_region) ;
    } else {
        load_def_region_files() ;
    }

    lprint "Loading Retime Related files : \n\n" ;
    load_port_map_file () ;
    load_routeRules_files () ;
    lprint "\nDone.\n\n" ;


    # to get the feedthr hash for every 2 partitions 
    if ($s_type eq 'noscan' || $s_type eq 'feflat') {
        get_all_par_thrs () ;
    }

    my $proj = $ENV{NV_PROJECT} ;

    $main::opt_ultra = 0;
    if ( defined $opt_ultra) {
        $main::opt_ultra = $opt_ultra;
    }

    my $inter_path_only = 0;
    if ( defined $opt_inter_only) {
        lprint "# Analysis inter-chiplet paths Only \n";
        $inter_path_only = 1;
    }

    %main::share_variables = ();
    ### Define per retime stage distance
    if ($s_type eq 'ipo' || $s_type eq 'flat') {
        $main::share_variables{"RETIME_DISTANCE"} = 1200 ;
    } else {
        $main::share_variables{"RETIME_DISTANCE"} = 1000 ;
    }
    if ( defined $retime_distance) {
        $main::share_variables{"RETIME_DISTANCE"} = $retime_distance;
    }

    ### Histogram Constants
    $main::share_variables{"DIST_HISTOGRAM_MIN"} = 700;
    if ( defined $histogram_min) {
        $main::share_variables{"DIST_HISTOGRAM_MIN"} = $hist_min;
    }

    $main::share_variables{"DIST_HISTOGRAM_MAX"} = 2000;
    if ( defined $histogram_max) {
        $main::share_variables{"DIST_HISTOGRAM_MAX"} = $histogram_max;
    }

    $main::share_variables{"DIST_HISTOGRAM_SCALE"} = 100;
    if ( defined $histogram_step) {
        $main::share_variables{"DIST_HISTOGRAM_SCALE"} = $histogram_step;
    }

    ######################
    ## violation filter ##
    ######################

    if (!defined $opt_no_clr) {
        clear_vios;
        load_vios $vio_file;
    }

    my @vios = ();

    if (($top eq "nv_top") || ($top eq "nvs_top") && ($inter_path_only == 1)) {
        @vios = all_path_vios(-filter => "is_inter_chiplet");
    } else {
        @vios = all_path_vios;
    }

    if (scalar (@skip_clock) && !scalar(@only_clock)) {
        print "INFO: User defined skip clocks @skip_clock\n\n";
        my $clk_regexp   = "" ;
        my @clk_patterns = () ;
        if (scalar @skip_clock > 1) {
            @clk_patterns = map (glob2regex_txt ($_), @skip_clock) ;
        } else {
            my @skip_clocks = split (" ", $skip_clock[0]) ;
            @clk_patterns = map (glob2regex_txt ($_), @skip_clocks) ;
        }
        $clk_regexp   = join ("|", @clk_patterns) ;
        @main::selected_vios = all_vios (-filter => "end_clk !~ /$clk_regexp/") ;

        lprint "Skip clocks : $clk_regexp" if (defined $opt_debug) ;

    } elsif (scalar (@only_clock) && !scalar (@skip_clock)){
        print "INFO: User defined only clocks @only_clock\n\n";
        my $clk_regexp   = "" ;
        my @clk_patterns = () ;
        if (scalar @only_clock > 1) {
            @clk_patterns = map (glob2regex_txt ($_), @only_clock) ;
        } else {
            my @only_clocks = split (" ", $only_clock[0]) ;
            @clk_patterns = map (glob2regex_txt ($_), @only_clocks) ;
        }
        $clk_regexp   = join ("|", @clk_patterns) ;
        @main::selected_vios = all_vios (-filter => "end_clk =~ /$clk_regexp/") ;

        lprint "Only clocks : $clk_regexp" if (defined $opt_debug) ;

    } elsif (scalar (@only_clock) && scalar (@skip_clock)) {
        error "DO not define skip_clock and only_clock at the same time which has conflict!\n\n";
        @main::selected_vios = ();
    } else {
        @main::selected_vios = @vios;
    }
 
    $timing_corner = get_timing_corner;
    $timing_corner =~ s/v/v_max_si/g;
    $main::mender_delay_per_dist_corner = 0;
    $main::mender_delay_per_dist_corner = get_yaml_corner_attr ($timing_corner, mender_delay_per_dist);
    lprint "\n# using $mender_delay_per_dist_corner ns per 1000um for $timing_corner delay estimate \n\n";

    ######################
    ## use multi-thread ##
    ######################

    my $curr_vio_count = scalar(@main::selected_vios);
    lprint "Total Violation Num : $curr_vio_count\n" ;
    my $thread_count = 16 ;
    if ($curr_vio_count < 200000) {
        $thread_count = 4 ;
    } elsif ($curr_vio_count < 400000) {
        $thread_count = 8 ;
    } else {
        $thread_count = 16 ;
    }
    if (defined $multi_thread_num) {
        $thread_count = $multi_thread_num;
    }
    my $count_per_thread = floor( $curr_vio_count / $thread_count );
    my $remainders = $curr_vio_count % $thread_count;
    my @vio_cnt_range;
    foreach my $thread (0..($thread_count - 1)) {
        if ($thread < $remainders) {
            push (@vio_cnt_range , "1:".($count_per_thread + 1).":$thread:$thread_count") ;
        } else {
            push (@vio_cnt_range , "1:$count_per_thread:$thread:$thread_count") ;
        }
    }

    my $start_date = `date`;
    chomp $start_date;

    lprint ("    Use multi-thread $thread_count \n\n");
    lprint ("    $start_date\n\n");
 
    %main::unit_hier_mapping        = genUnitNameHash ;
    %main::p_feed_pars              = () ;
    %main::p_feed_pars_num          = () ;
    %main::p_start_unit             = () ;
    %main::p_end_unit               = () ;
    %main::p_par_num                = () ;
    %main::p_module_split_by_par    = () ;
    %main::p_man_dist               = () ;
    %main::p_real_dist              = () ;
    %main::p_ideal_dist             = () ;
    %main::p_mcp_setup_num          = () ;
    %main::p_detour_ratio           = () ;
    %main::p_bin_man_dist           = () ;
    %main::p_bin_ideal_dist         = () ;
    %main::p_bin_real_dist          = () ;
    %main::p_is_mcp                 = () ;
    %main::p_sig_name               = () ;
    %main::p_end_clk                = () ;
    %main::p_clk_period             = () ;
    %main::p_start_routeRule        = () ;
    %main::p_end_routeRule          = () ;
    %main::p_source_coor            = () ;
    %main::p_dest_coor              = () ;
    %main::p_source_port_dist       = () ;
    %main::p_dest_port_dist         = () ;

    # info for dumping reps
    %main::p_dist_rep               = () ;
    %main::p_io_dist_rep            = () ;
    %main::p_mann_plot_rep          = () ;
    %main::p_detour_plot_rep        = () ;
    %main::p_retime_rep             = () ;
    %main::p_mcp_rep                = () ;
    %main::p_rd_sum_rep             = () ;

    %main::p_AdjParNRT2NRT_mcp      = () ;
    %main::p_AdjParNRT2NRT_nonmcp   = () ;
    %main::p_AdjParRT2NRT_mcp       = () ;
    %main::p_AdjParRT2NRT_nonmcp    = () ;
    %main::p_AdjParNRT2RT_mcp       = () ;
    %main::p_AdjParNRT2RT_nonmcp    = () ;
    %main::p_AdjParRT2RT_mcp        = () ;
    %main::p_AdjParRT2RT_nonmcp     = () ;

    %main::p_max_feed_pars_num      = () ;
    %main::p_max_feed_pars          = () ;
    %main::p_max_man_dist           = () ;
    %main::p_max_ideal_dist         = () ;
    %main::p_max_real_dist          = () ;
    %main::p_max_source_port_dist   = () ;
    %main::p_max_dest_port_dist     = () ;
    %main::p_max_end_clk            = () ;
    %main::p_sig_is_mcp             = () ;
    %main::p_sig_s_routeRule        = () ;
    %main::p_sig_e_routeRule        = () ;

    %main::p_max_nonrt_mcp_feed_pars_num    = () ;
    %main::p_max_nonrt_mcp_feed_pars        = () ;
    %main::p_max_nonrt_mcp_man_dist         = () ;
    %main::p_max_nonrt_mcp_ideal_dist       = () ;
    %main::p_max_nonrt_mcp_real_dist        = () ;
    %main::p_max_nonrt_mcp_source_port_dist = () ;
    %main::p_max_nonrt_mcp_dest_port_dist   = () ;
    %main::p_max_nonrt_mcp_end_clk          = () ;
    %main::p_max_nonrt_mcp_clks             = () ;
    %main::p_max_nonrt_mcp_worst_id         = () ;
    %main::p_nonrt_mcp_sig_is_mcp           = () ;

    %main::p_max_nonrt_nonmcp_feed_pars_num    = () ;
    %main::p_max_nonrt_nonmcp_feed_pars        = () ;
    %main::p_max_nonrt_nonmcp_man_dist         = () ;
    %main::p_max_nonrt_nonmcp_ideal_dist       = () ;
    %main::p_max_nonrt_nonmcp_real_dist        = () ;
    %main::p_max_nonrt_nonmcp_source_port_dist = () ;
    %main::p_max_nonrt_nonmcp_dest_port_dist   = () ;
    %main::p_max_nonrt_nonmcp_end_clk          = () ;
    %main::p_max_nonrt_nonmcp_clks             = () ;
    %main::p_max_nonrt_nonmcp_worst_id         = () ;
    %main::p_nonrt_nonmcp_sig_is_mcp           = () ;

    %main::p_max_rt_mcp_feed_pars_num    = () ;
    %main::p_max_rt_mcp_feed_pars        = () ;
    %main::p_max_rt_mcp_man_dist         = () ;
    %main::p_max_rt_mcp_ideal_dist       = () ;
    %main::p_max_rt_mcp_real_dist        = () ;
    %main::p_max_rt_mcp_source_port_dist = () ;
    %main::p_max_rt_mcp_dest_port_dist   = () ;
    %main::p_max_rt_mcp_end_clk          = () ;
    %main::p_max_rt_mcp_s_routeRule      = () ;
    %main::p_max_rt_mcp_e_routeRule      = () ;
    %main::p_rt_mcp_sig_is_mcp           = () ;
    %main::p_max_rt_mcp_clks             = () ;
    %main::p_max_rt_mcp_worst_id         = () ;

    %main::p_max_rt_nonmcp_feed_pars_num    = () ;
    %main::p_max_rt_nonmcp_feed_pars        = () ;
    %main::p_max_rt_nonmcp_man_dist         = () ;
    %main::p_max_rt_nonmcp_ideal_dist       = () ;
    %main::p_max_rt_nonmcp_real_dist        = () ;
    %main::p_max_rt_nonmcp_source_port_dist = () ;
    %main::p_max_rt_nonmcp_dest_port_dist   = () ;
    %main::p_max_rt_nonmcp_end_clk          = () ;
    %main::p_max_rt_nonmcp_s_routeRule      = () ;
    %main::p_max_rt_nonmcp_e_routeRule      = () ;
    %main::p_rt_nonmcp_sig_is_mcp           = () ;
    %main::p_max_rt_nonmcp_clks             = () ;
    %main::p_max_rt_nonmcp_worst_id         = () ;


    %main::p_FeedParNRT2NRT_mcp     = () ;
    %main::p_FeedParNRT2NRT_nonmcp  = () ;
    %main::p_FeedParRT_mcp          = () ;
    %main::p_FeedParRT_nonmcp       = () ;

    %main::p_dist_histogram_mann    = () ;
    %main::p_dist_histogram_ideal   = () ;
    %main::p_dist_histogram_real    = () ;
    %main::p_detour_histogram       = () ;

    %main::p_comb_detour_start_routeRule  = () ;
    %main::p_comb_detour_end_routeRule    = () ;
    %main::p_comb_detour_feed_pars_num    = () ;
    %main::p_comb_detour_feed_pars        = () ;
    %main::p_comb_detour_man_dist         = () ;
    %main::p_comb_detour_ideal_dist       = () ;
    %main::p_comb_detour_real_dist        = () ;
    %main::p_comb_detour_source_port_dist = () ;
    %main::p_comb_detour_dest_port_dist   = () ;
    %main::p_comb_detour_end_clk          = () ;
    %main::p_comb_detour_is_mcp           = () ;

    %main::p_intra_par_vio                = () ;

    %main::p_io_vio                       = () ;
    %main::p_io_max_vio                   = () ;
    %main::p_io_worst_id                  = () ;

    @main::all_partitions      = () ;
    if ($top ne "nv_top" && $top ne "nvs_top") {
        @main::all_partitions = get_all_par_insts();
    } else {
         @main::all_partitions = ();
    }
    
    alarm (0); #Avoid heartbeat interference with fork manager
    my $pl = Parallel::Loops->new($thread_count);
    
    # share varibles/array/hash in multi-threads
    $pl->share(\%main::unit_hier_mapping);
    $pl->share(\%main::p_feed_pars);
    $pl->share(\%main::p_feed_pars_num);
    $pl->share(\%main::p_start_unit);
    $pl->share(\%main::p_end_unit);
    $pl->share(\%main::p_par_num);
    $pl->share(\%main::p_module_split_by_par);
    $pl->share(\%main::p_man_dist);
    $pl->share(\%main::p_real_dist);
    $pl->share(\%main::p_ideal_dist);
    $pl->share(\%main::p_detour_ratio);
    $pl->share(\%main::p_mcp_setup_num);
    $pl->share(\%main::p_bin_man_distance);
    $pl->share(\%main::p_bin_ideal_distance);
    $pl->share(\%main::p_bin_real_distance);
    $pl->share(\@main::selected_vios);
    $pl->share(\%main::share_variables);
    $pl->share(\%main::p_is_mcp);
    $pl->share(\%main::p_sig_name) ;
    $pl->share(\%main::p_end_clk) ;
    $pl->share(\%main::p_clk_period) ;
    $pl->share(\%main::p_start_routeRule) ;
    $pl->share(\%main::p_end_routeRule) ;
    $pl->share(\%main::p_source_coor) ;
    $pl->share(\%main::p_dest_coor) ;
    $pl->share(\%main::p_source_port_dist) ;
    $pl->share(\%main::p_dest_port_dist) ;
    $pl->share(\%main::p_dist_rep) ;
    $pl->share(\%main::p_io_dist_rep) ;
    $pl->share(\%main::p_dist_rep) ;
    $pl->share(\%main::p_io_dist_rep) ;
    $pl->share(\%main::p_mann_plot_rep) ;
    $pl->share(\%main::p_detour_plot_rep) ;
    $pl->share(\%main::p_retime_rep) ;
    $pl->share(\%main::p_mcp_rep) ;
    $pl->share(\%main::p_rd_sum_rep) ;

    $pl->share(\%main::p_AdjParNRT2NRT_mcp) ;
    $pl->share(\%main::p_AdjParNRT2NRT_nonmcp) ;
    $pl->share(\%main::p_AdjParRT2NRT_mcp) ;
    $pl->share(\%main::p_AdjParRT2NRT_nonmcp) ;
    $pl->share(\%main::p_AdjParNRT2RT_mcp) ;
    $pl->share(\%main::p_AdjParNRT2RT_nonmcp) ;
    $pl->share(\%main::p_AdjParRT2RT_mcp) ;
    $pl->share(\%main::p_AdjParRT2RT_nonmcp) ;

    $pl->share(\%main::p_max_feed_pars_num) ;
    $pl->share(\%main::p_max_feed_pars) ;
    $pl->share(\%main::p_max_feed_pars) ;
    $pl->share(\%main::p_max_man_dist) ;
    $pl->share(\%main::p_max_ideal_dist) ;
    $pl->share(\%main::p_max_real_dist) ;
    $pl->share(\%main::p_max_source_port_dist) ;
    $pl->share(\%main::p_max_dest_port_dist) ;
    $pl->share(\%main::p_max_end_clk) ;
    $pl->share(\%main::p_sig_is_mcp) ;
    $pl->share(\%main::p_sig_s_routeRule) ;
    $pl->share(\%main::p_sig_e_routeRule) ;
    $pl->share(\%main::p_FeedParNRT2NRT_mcp) ;
    $pl->share(\%main::p_FeedParNRT2NRT_nonmcp) ;
    $pl->share(\%main::p_FeedParRT_mcp) ;
    $pl->share(\%main::p_FeedParRT_nonmcp) ;

    $pl->share(\%main::p_max_nonrt_mcp_feed_pars_num) ;
    $pl->share(\%main::p_max_nonrt_mcp_feed_pars) ;
    $pl->share(\%main::p_max_nonrt_mcp_man_dist) ;
    $pl->share(\%main::p_max_nonrt_mcp_ideal_dist) ;
    $pl->share(\%main::p_max_nonrt_mcp_real_dist) ;
    $pl->share(\%main::p_max_nonrt_mcp_source_port_dist) ;
    $pl->share(\%main::p_max_nonrt_mcp_dest_port_dist) ;
    $pl->share(\%main::p_max_nonrt_mcp_end_clk) ;
    $pl->share(\%main::p_nonrt_mcp_sig_is_mcp) ;
    $pl->share(\%main::p_max_nonrt_mcp_clks) ;
    $pl->share(\%main::p_max_nonrt_mcp_worst_id) ;

    $pl->share(\%main::p_max_nonrt_nonmcp_feed_pars_num) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_feed_pars) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_man_dist) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_ideal_dist) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_real_dist) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_source_port_dist) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_dest_port_dist) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_end_clk) ;
    $pl->share(\%main::p_nonrt_nonmcp_sig_is_mcp) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_clks) ;
    $pl->share(\%main::p_max_nonrt_nonmcp_worst_id) ;

    $pl->share(\%main::p_max_rt_mcp_feed_pars_num) ;
    $pl->share(\%main::p_max_rt_mcp_feed_pars) ;
    $pl->share(\%main::p_max_rt_mcp_man_dist) ;
    $pl->share(\%main::p_max_rt_mcp_ideal_dist) ;
    $pl->share(\%main::p_max_rt_mcp_real_dist) ;
    $pl->share(\%main::p_max_rt_mcp_source_port_dist) ;
    $pl->share(\%main::p_max_rt_mcp_dest_port_dist) ;
    $pl->share(\%main::p_max_rt_mcp_end_clk) ;
    $pl->share(\%main::p_max_rt_mcp_s_routeRule) ;
    $pl->share(\%main::p_max_rt_mcp_e_routeRule) ;
    $pl->share(\%main::p_rt_mcp_sig_is_mcp) ;
    $pl->share(\%main::p_max_rt_mcp_clks) ;
    $pl->share(\%main::p_max_rt_mcp_worst_id) ;

    $pl->share(\%main::p_max_rt_nonmcp_feed_pars_num) ;
    $pl->share(\%main::p_max_rt_nonmcp_feed_pars) ;
    $pl->share(\%main::p_max_rt_nonmcp_man_dist) ;
    $pl->share(\%main::p_max_rt_nonmcp_ideal_dist) ;
    $pl->share(\%main::p_max_rt_nonmcp_real_dist) ;
    $pl->share(\%main::p_max_rt_nonmcp_source_port_dist) ;
    $pl->share(\%main::p_max_rt_nonmcp_dest_port_dist) ;
    $pl->share(\%main::p_max_rt_nonmcp_end_clk) ;
    $pl->share(\%main::p_max_rt_nonmcp_s_routeRule) ;
    $pl->share(\%main::p_max_rt_nonmcp_e_routeRule) ;
    $pl->share(\%main::p_rt_nonmcp_sig_is_mcp) ;
    $pl->share(\%main::p_max_rt_nonmcp_clks) ;
    $pl->share(\%main::p_max_rt_nonmcp_worst_id) ;

    $pl->share(\%main::p_dist_histogram_mann) ;
    $pl->share(\%main::p_dist_histogram_ideal) ;
    $pl->share(\%main::p_dist_histogram_real) ;
    $pl->share(\%main::p_detour_histogram) ;

    $pl->share(\%main::p_comb_detour_start_routeRule) ;
    $pl->share(\%main::p_comb_detour_end_routeRule) ;
    $pl->share(\%main::p_comb_detour_feed_pars_num) ;
    $pl->share(\%main::p_comb_detour_feed_pars) ;
    $pl->share(\%main::p_comb_detour_man_dist) ;
    $pl->share(\%main::p_comb_detour_ideal_dist) ;
    $pl->share(\%main::p_comb_detour_real_dist) ;
    $pl->share(\%main::p_comb_detour_source_port_dist) ;
    $pl->share(\%main::p_comb_detour_dest_port_dist) ;
    $pl->share(\%main::p_comb_detour_end_clk) ;
    $pl->share(\%main::p_comb_detour_is_mcp) ;
    $pl->share(\%main::p_comb_detour_clks) ;
    $pl->share(\%main::p_comb_detour_worst_id) ;

    $pl->share(\%main::p_intra_par_vio) ;

    $pl->share(\%main::p_io_vio) ;
    $pl->share(\%main::p_io_max_vio) ;
    $pl->share(\%main::p_io_worst_id) ;

    #$s_type = session_type() ;
    #chomp $s_type ;

    ## Initialize Histogram
    my $loopCnt;
    for ($loopCnt = $main::share_variables{DIST_HISTOGRAM_MIN} ; $loopCnt < $main::share_variables{DIST_HISTOGRAM_MAX} ; $loopCnt += $main::share_variables{DIST_HISTOGRAM_SCALE}) {
        $main::p_dist_histogram_mann{$loopCnt}  = 0 ;
        $main::p_dist_histogram_ideal{$loopCnt} = 0 ;
        $main::p_dist_histogram_real{$loopCnt}  = 0 ;
    }
    $main::p_dist_histogram_mann{$main::share_variables{DIST_HISTOGRAM_MAX}}  = 0 ; 
    $main::p_dist_histogram_ideal{$main::share_variables{DIST_HISTOGRAM_MAX}} = 0 ; 
    $main::p_dist_histogram_real{$main::share_variables{DIST_HISTOGRAM_MAX}}  = 0 ; 

    for ($loopCnt = 1.0 ; $loopCnt <= 2.0 ; $loopCnt += 0.1) {
        $main::p_detour_histogram{$loopCnt} = 0 ;
    }
    $main::p_detour_histogram{2.0} = 0 ;
    

    # Kickoff multi-thread runs
    $pl->foreach(\@vio_cnt_range, \&cal_vio_attribute_by_index);
 
    my $end_date = `date`;
    chomp $end_date;
    lprint ("\n    Start Cal Attribute $start_date\n");
    lprint ("    Start Tag Attribute $end_date\n");

    if (defined $opt_dump_mcp_tcl) {
        dump_mcp_tcl_file ($rep_name) ;
    }

    dump_reports($top,$rep_name);

    my $subject = $rep_name ;
    $subject =~ s/\S+\/(\S+)\/\S+/$1/ ;
    $subject = "Retime/Detour Check $subject done." ;
    my $body = $rep_name ;
    $body =~ s/(\S+\/)\S+/$1/ ; 
    $body = "Reports @ :\n\t$body" ;

    if (defined $mail_id) {
        my $cmd  = "echo \"$body\" | mutt $mail_id -s \"$subject\"" ;
        system "$cmd" ;
        lprint "$body" ; 
    }
    

    undef $debug ;
    undef $log_file ;
    close_log ;

END

sub cal_vio_attribute_by_index {
    my ($index) = $_;
    my ($star_index,$end_index,$curr_thread,$thread_count) = split (/:/, $index);
    my $OUT = "OUT${curr_thread}" ;
    if ($debug) {
        my $debug_log = $log_file ;
        $debug_log =~ s/\.log$/_thread${curr_thread}.log/ ;
        open "$OUT" , "> $debug_log" or die "Can't write to $debug_log\n" ;
    }
    my $max_id = scalar(@main::selected_vios);
    foreach my $id ($star_index..$end_index) {
        my $real_index = ($id - 1)*$thread_count + $curr_thread;
        my $status_msg ;
        if ($id > 10000) {
            $status_msg = ($id) % 10000 ;
        } else {
            $status_msg = ($id) % 1000 ;
        }
        if ($status_msg == 0) {
            lprint("    Processed  $id of $end_index in Thread $curr_thread  \n");
        }
        if ($id == $end_index) {
            lprint("    Processed All $end_index in Thread $curr_thread  \n");
        }
        if (defined $debug) {
            &cal_vio_attribute($main::selected_vios[$real_index], $OUT);
        } else {
            &cal_vio_attribute($main::selected_vios[$real_index]) ;
        }
    }
    if ($debug) {
        close "$OUT" ;
    }
}

sub cal_vio_attribute {
    my $curr_vio ;
    my $OUT ;

    if ($debug) {
        ($curr_vio, $OUT) = @_ ;
    } else {
        ($curr_vio) = @_ ;
    }
 

    my $id         = attr_of_path_vio(id,$curr_vio);
    my $top        = get_top;
    my $startpoint = attr_of_path_vio(start_pin, $curr_vio);
    $startpoint    =~ s/checkpin.*//g; # this line is in order to replace the checkin(internal pin for analog cell such as PLL) due to mender can't see it 
    my $endpoint   = attr_of_path_vio(end_pin, $curr_vio);


    my @feed_pars            = ();
    my $capt_clk             = attr_of_vio ('end_clk' => $curr_vio);
    my $period               = attr_of_vio ('period' => $curr_vio);
    my $capture_time         = attr_of_vio ('end_clk_edge_dly' => $curr_vio);
    my $endpar               = attr_of_vio ('end_par' => $curr_vio);
    my $startpar             = attr_of_vio ('start_par' => $curr_vio);
    my $is_io                = attr_of_vio ('is_io' => $curr_vio);
    my @module_split_by_par  = attr_of_vio ('module_split_by_par' => $curr_vio);
    my $module_split_by_par  = join ('->', @module_split_by_par);
    my $num                  = scalar @module_split_by_par;
    my $start_unit           = attr_of_vio ('start_unit' => $curr_vio) ;
    my $end_unit             = attr_of_vio ('end_unit' => $curr_vio) ;
    my @vio_thr_pars         = get_vio_thr_pars ($curr_vio) ;


    $main::p_start_unit{$id} = $start_unit ; 
    $main::p_end_unit{$id}   = $end_unit ; 
    $main::p_par_num{$id}    = $num;
    $main::p_sig_name{$id}   = GetVioSigName($curr_vio) ;
    $main::p_end_clk{$id}    = $capt_clk ;
    $main::p_clk_period{$id} = $period ; 
 
    $main::p_module_split_by_par{$id} = $module_split_by_par;

    if ($debug) {
        print $OUT "VIO : $id $startpoint $endpoint $capt_clk $main::p_sig_name{$id}\n" ;
    }

    # global variable $top and $s_type
    # $s_type = session_type() ;
    # chomp $s_type ;

    if ( ($top ne "nv_top") && ($top ne "nvs_top") && ($s_type eq "noscan" | $s_type eq "feflat")) {
        if($main::p_sig_name{$id} ne "NA") {
            @feed_pars = get_vio_feeds($curr_vio) ;
            $feed_pars_num = scalar @feed_pars;
            $feed_pars = join('->',@feed_pars);
            $main::p_feed_pars{$id} = $feed_pars;
            $main::p_feed_pars_num{$id} = $feed_pars_num;
        } else {
            $main::p_feed_pars{$id} = join('->',@vio_thr_pars);
            $main::p_feed_pars_num{$id} = $num;
        }
    } else {
        $main::p_feed_pars{$id} = join('->',@vio_thr_pars);
        $main::p_feed_pars_num{$id} = $num;
    }

    ### lat/D pin as startpoint will fail real_dist()
    $startpoint_lat = attr_of_pin(is_latch, $startpoint);
    if ($startpoint_lat == 1) {
        $startpoint =~ s/\/D$/\/Q/;
    }
    my $mann_dist = get_dist ($startpoint => $endpoint);
    my @pin_list = split (/ /, attr_of_path_vio(pin_list, $curr_vio));
    map ($_ =~ s/checkpin.*//g, @pin_list); # this line is in order to replace the checkin(internal pin for analog cell such as PLL) due to mender can't see it
    my $real_dist = get_dist_pinArray (@pin_list);
    @pin_list = $startpoint;
    push (@pin_list => attr_of_path_vio(inter_pin_array, $curr_vio));
    push (@pin_list => $endpoint);
    my $ideal_dist = get_dist_pinArray (@pin_list);
    my $detour = -999;
    if ( $mann_dist > 0 ) {
        $detour = $ideal_dist / $mann_dist ;
    }
    my $detour_ratio = sprintf('%.2f',$detour);

    $main::p_detour{$id}       = $detour;
    $main::p_detour_ratio{$id} = $detour_ratio;
    $main::p_man_dist{$id}     = $mann_dist;
    $main::p_real_dist{$id}    = $real_dist;
    $main::p_ideal_dist{$id}   = $ideal_dist;

    my $mcp_num = 0;
    my $real_setup_mcp_num = 1;
    my $real_hold_mcp_num  = 0;

    if ($main::opt_ultra == 1) {
        $mcp_num = ($ideal_dist / 1000) * $main::mender_delay_per_dist_corner / $period;
        $real_setup_mcp_num = ceil($capture_time / $period);
        $real_hold_mcp_num = $real_setup_mcp_num - 1;
    } else {
        $mcp_num = 0;
        #fix precision issue in clock period and capture_time
        $capture_time = int($capture_time*100);
        $period = int($period*100);
        if ($capture_time > $period) {
            $real_setup_mcp_num = 2;
        }
        $real_hold_mcp_num = 1;
    }

    my $mcp_hold_num  = int $mcp_num;
    my $mcp_setup_num = $mcp_hold_num + 1;

    $main::p_mcp_setup_num{$id} = $mcp_setup_num;
    $main::p_mcp_num{$id}       =  $mcp_num;
    if ($real_setup_mcp_num > 1 ) {
        $main::p_is_mcp{$id} = 1;
    } else {
        $main::p_is_mcp{$id} = 0;
    }

    my $bin_man_dist   = (int ($man_dist / 100) + 1) * 100;
    my $bin_ideal_dist = (int ($ideal_dist / 100) + 1) * 100;
    my $bin_real_dist  = (int ($real_dist / 100) + 1) * 100;
    $main::p_bin_man_distance{$id}   = $bin_man_dist;
    $main::p_bin_ideal_distance{$id} = $bin_ideal_dist;
    $main::p_bin_real_distance{$id}  = $bin_real_dist;
    
    #$main::p_start_routeRule{$id}  = GetStartRouteRule($curr_vio) ;
    #$main::p_end_routeRule{$id}    = GetEndRouteRule($curr_vio) ;
    $main::p_start_routeRule{$id}  = "NA" ;
    $main::p_end_routeRule{$id}    = "NA" ;

    # the info for dumping reports

    $main::p_dist_rep{$id} = "$startpoint -> $endpoint, Manhattan: $mann_dist, Ideal: $ideal_dist, Real: $real_dist, detour: $main::p_detour_ratio{$id}, is_mcp: $main::p_is_mcp{$id} , clk: $capt_clk";

    if ($is_io) {
        $main::p_io_dist_rep{$id} = "$startpoint -> $endpoint, Manhattan: $mann_dist, Ideal: $ideal_dist, Real: $real_dist, detour: $main::p_detour_ratio{$id}, clk: $capt_clk";
    }

    if ( $mann_dist >= $main::share_variables{"DIST_HISTOGRAM_MIN"} ) {
        $main::p_mann_plot_rep{$id} = "\@tmp1 = get_pin_xy ${startpoint}; \@tmp2 = get_pin_xy ${endpoint} ; plot_line -color red \@tmp1 \@tmp2; ## Manhattan: $mann_dist, clk: $capt_clk," ;
    }
 
    if ( $ideal_dist >= $main::share_variables{RETIME_DISTANCE}) {
        # need check chiplet level real_distance > RETIME_DISTANCE and feed_pars >=3 paths
        # need check nv_top(since all paths are inter-chiplet), man_distance > RETIME_DISTANCE
        if ($main::p_feed_pars_num{$id} >= 3 && $main::p_sig_name{$id} ne "" && $s_type ne 'ipo' && $s_type ne 'flat') {
            if ($mann_dist <= $main::share_variables{RETIME_DISTANCE} && ($s_type eq "noscan" | $s_type eq "feflat")) {
                $main::p_retime_rep{$id} = "ATTENTION : $main::p_sig_name{$id}\t\t$capt_clk\t$main::p_feed_pars{$id}\t\t$mann_dist\t$ideal_dist\tis_mcp: $main::p_is_mcp{$id}" ;
            } else {
                $main::p_retime_rep{$id} = "$main::p_sig_name{$id}\t\t$capt_clk\t$main::p_feed_pars{$id}\t\t$mann_dist\t$ideal_dist\tis_mcp: $main::p_is_mcp{$id}" ;
            }
        } elsif ($main::p_feed_pars_num{$id} >= 2 && $main::p_sig_name{$id} ne "" && ($s_type eq 'ipo' || $s_type eq 'flat')) {
            $main::p_retime_rep{$id} = "$main::p_sig_name{$id}\t\t$capt_clk\t$main::p_feed_pars{$id}\t\t$mann_dist\t$ideal_dist\tis_mcp: $main::p_is_mcp{$id}" ;
        } elsif (($top eq "nv_top") || ($top eq "nvs_top")) {
        # deal with those feed_par_num == 2's inter-chiplet paths(noscan, if there's no combination, feedpars is 2 for most inter-chiplet paths)
            if ($mann_dist >= $main::share_variables{RETIME_DISTANCE}*1.5) {
                $main::p_retime_rep{$id} = "$main::p_sig_name{$id}\t\t$capt_clk\t$main::p_feed_pars{$id}\t\t$mann_dist\t$ideal_dist\tis_mcp: $main::p_is_mcp{$id}" ;
            }
        }
    }

    if ( $ideal_dist >= $main::share_variables{RETIME_DISTANCE}) {
        if ($mcp_hold_num > 0 && $main::p_feed_pars_num{$id} >= 3) {
            if ($real_setup_mcp_num > 1 ) {
                $main::p_mcp_rep{$id} = "# real mcp is $real_setup_mcp_num: nv_set_mcp -setup $mcp_setup_num -from [get_pins $startpoint] -to [get_pins $endpoint] -infor {NV_ERR: this is temporary for retime missing} \n" ;
                $main::p_mcp_rep{$id} = $main::p_mcp_rep{$id} . "# real mcp is $real_hold_mcp_num: nv_set_mcp -hold $mcp_hold_num -from [get_pins $startpoint] -to [get_pins $endpoint] -infor {NV_ERR: this is temporary for retime missing}" ;
            } else {
                $main::p_mcp_rep{$id} = "nv_set_mcp -setup $mcp_setup_num -from [get_pins $startpoint] -to [get_pins $endpoint] -infor {NV_ERR: this is temporary for retime missing} \n" ;
                $main::p_mcp_rep{$id} = $main::p_mcp_rep{$id} . "nv_set_mcp -hold $mcp_hold_num -from [get_pins $startpoint] -to [get_pins $endpoint] -infor {NV_ERR: this is temporary for retime missing} " ;
            }
        }
    }

    if ( $main::p_detour_ratio{$id} >= 1.2 && $main::p_ideal_dist{$id} >= $main::share_variables{"DIST_HISTOGRAM_MIN"} ) {
        $outStr = join(' -through ', @pin_list);
        $main::p_detour_plot_rep{$id} = "mrt \-plot -from ${outStr} ; ## detour: $main::p_detour_ratio{$id}, Ideal: $ideal_dist, clk: $capt_clk,";
    }

    if ($debug) {
        print $OUT "SIG_NAME : $startpar $endpar $main::p_sig_name{$id} \n";
    }

    # to fix, not startpar, but start_par_cell 
    my $start_par_cell = attr_of_vio ('start_par_cell' => $curr_vio) ;
    my $end_par_cell   = attr_of_vio ('end_par_cell' => $curr_vio) ;
    if ($start_par_cell eq $end_par_cell && $main::p_sig_name{$id} ne "NA" && $main::p_feed_pars_num{$id} > 1 && $is_io == 0 && $main::p_is_mcp{$id} == 0) {
        $main::p_start_routeRule{$id}  = GetStartRouteRule($curr_vio) ;
        $main::p_end_routeRule{$id}    = GetEndRouteRule($curr_vio) ;

        $main::p_comb_detour_start_routeRule{$main::p_sig_name{$id}}  = $main::p_start_routeRule{$id} ;
        $main::p_comb_detour_end_routeRule{$main::p_sig_name{$id}}    = $main::p_end_routeRule{$id} ;
        $main::p_comb_detour_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
        $main::p_comb_detour_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
        $main::p_comb_detour_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
        $main::p_comb_detour_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
        $main::p_comb_detour_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
        $main::p_comb_detour_end_clk{$main::p_sig_name{$id}}          = $capt_clk ;
        $main::p_comb_detour_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
        $main::p_comb_detour_worst_id{$main::p_sig_name{$id}}         = $id ;
    }

    my $dist_threshold = $main::share_variables{RETIME_DISTANCE} ;

    # enlarge the distance threshold for slow clocks
    if ($main::p_clk_period{$id} >= 2) {
        $dist_threshold = 750 * $main::p_clk_period{$id} ;
    }
    if ($dist_threshold < $main::share_variables{RETIME_DISTANCE}) {
        $dist_threshold = $main::share_variables{RETIME_DISTANCE} ;
    }

    # the chiplet port target should be half of the one inside chiplet
    if ($is_io) {
        $dist_threshold = $dist_threshold * 0.5 ;
    }
    
  
    if ($ideal_dist >= $dist_threshold) {
        # feedthr pars ------------------- nonRT => nonRT # missing retiming
        #                    |------------- RT related    # error 
        if ($main::p_feed_pars_num{$id} > 2 && $main::p_sig_name{$id} ne "NA") {
            if ($start_unit !~ /_retime_partition_/ && $end_unit !~ /_retime_partition_/) {
                if ($main::p_is_mcp{$id}) {
                    if (exists $main::p_max_nonrt_mcp_feed_pars_num{$main::p_sig_name{$id}}) {
                        if ($main::p_feed_pars_num{$id} > $main::p_max_nonrt_mcp_feed_pars_num{$main::p_sig_name{$id}}) {
                            $main::p_max_nonrt_mcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                            $main::p_max_nonrt_mcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                            $main::p_max_nonrt_mcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                            $main::p_max_nonrt_mcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                            $main::p_max_nonrt_mcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                            $main::p_max_nonrt_mcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                            $main::p_nonrt_mcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                            $main::p_max_nonrt_mcp_worst_id{$main::p_sig_name{$id}} = $id ; 
                        }
                    } else {
                        $main::p_max_nonrt_mcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                        $main::p_max_nonrt_mcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                        $main::p_max_nonrt_mcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                        $main::p_max_nonrt_mcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                        $main::p_max_nonrt_mcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                        $main::p_max_nonrt_mcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                        $main::p_nonrt_mcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                        $main::p_max_nonrt_mcp_worst_id{$main::p_sig_name{$id}} = $id ; 
                    }
                } else {
                    if (exists $main::p_max_nonrt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}) {
                        if ($main::p_feed_pars_num{$id} > $main::p_max_nonrt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}) {
                            $main::p_max_nonrt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                            $main::p_max_nonrt_nonmcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                            $main::p_max_nonrt_nonmcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                            $main::p_max_nonrt_nonmcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                            $main::p_max_nonrt_nonmcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                            $main::p_max_nonrt_nonmcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                            $main::p_max_nonrt_nonmcp_worst_id{$main::p_sig_name{$id}} = $id ; 
                            $main::p_nonrt_nonmcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                        }
                    } else { 
                        $main::p_max_nonrt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                        $main::p_max_nonrt_nonmcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                        $main::p_max_nonrt_nonmcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                        $main::p_max_nonrt_nonmcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                        $main::p_max_nonrt_nonmcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                        $main::p_max_nonrt_nonmcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                        $main::p_max_nonrt_nonmcp_worst_id{$main::p_sig_name{$id}} = $id ; 
                        $main::p_nonrt_nonmcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                    }        
                }
            } else {
                $main::p_start_routeRule{$id}  = GetStartRouteRule($curr_vio) ;
                $main::p_end_routeRule{$id}    = GetEndRouteRule($curr_vio) ;
                if ($main::p_is_mcp{$id}) {
                    if (exists $main::p_max_rt_mcp_feed_pars_num{$main::p_sig_name{$id}}) {
                        if ($main::p_feed_pars_num{$id} > $main::p_max_rt_mcp_feed_pars_num{$main::p_sig_name{$id}}) {
                            $main::p_max_rt_mcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                            $main::p_max_rt_mcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                            $main::p_max_rt_mcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                            $main::p_max_rt_mcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                            $main::p_max_rt_mcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                            $main::p_max_rt_mcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                            $main::p_max_rt_mcp_s_routeRule{$main::p_sig_name{$id}}      = $main::p_start_routeRule{$id} ; 
                            $main::p_max_rt_mcp_e_routeRule{$main::p_sig_name{$id}}      = $main::p_end_routeRule{$id} ; 
                            $main::p_rt_mcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                            $main::p_max_rt_mcp_worst_id{$main::p_sig_name{$id}}         = $id ; 
                        }
                    } else { 
                        $main::p_max_rt_mcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                        $main::p_max_rt_mcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                        $main::p_max_rt_mcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                        $main::p_max_rt_mcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                        $main::p_max_rt_mcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                        $main::p_max_rt_mcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                        $main::p_max_rt_mcp_s_routeRule{$main::p_sig_name{$id}}      = $main::p_start_routeRule{$id} ; 
                        $main::p_max_rt_mcp_e_routeRule{$main::p_sig_name{$id}}      = $main::p_end_routeRule{$id} ; 
                        $main::p_rt_mcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                        $main::p_max_rt_mcp_worst_id{$main::p_sig_name{$id}}         = $id ; 
                    }     
                } else {
                    if (exists $main::p_max_rt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}) {
                        if ($main::p_feed_pars_num{$id} > $main::p_max_rt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}) {
                            $main::p_max_rt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                            $main::p_max_rt_nonmcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                            $main::p_max_rt_nonmcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                            $main::p_max_rt_nonmcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                            $main::p_max_rt_nonmcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                            $main::p_max_rt_nonmcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                            $main::p_max_rt_nonmcp_s_routeRule{$main::p_sig_name{$id}}      = $main::p_start_routeRule{$id} ; 
                            $main::p_max_rt_nonmcp_e_routeRule{$main::p_sig_name{$id}}      = $main::p_end_routeRule{$id} ; 
                            $main::p_rt_nonmcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                            $main::p_max_rt_nonmcp_worst_id{$main::p_sig_name{$id}}         = $id ; 
                        }
                    } else { 
                        $main::p_max_rt_nonmcp_feed_pars_num{$main::p_sig_name{$id}}    = $main::p_feed_pars_num{$id} ;
                        $main::p_max_rt_nonmcp_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
                        $main::p_max_rt_nonmcp_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
                        $main::p_max_rt_nonmcp_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
                        $main::p_max_rt_nonmcp_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
                        $main::p_max_rt_nonmcp_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
                        $main::p_max_rt_nonmcp_s_routeRule{$main::p_sig_name{$id}}      = $main::p_start_routeRule{$id} ; 
                        $main::p_max_rt_nonmcp_e_routeRule{$main::p_sig_name{$id}}      = $main::p_end_routeRule{$id} ; 
                        $main::p_rt_nonmcp_sig_is_mcp{$main::p_sig_name{$id}}           = $main::p_is_mcp{$id} ;
                        $main::p_max_rt_nonmcp_worst_id{$main::p_sig_name{$id}}         = $id ; 
                    }     
                }
            }
        # adjacent pars (ipo or flat)  ------ nonRT => nonRT # missing retiming 
        #                              |----- nonRT => RT    # start_port_dist/end_port_dist too large
        #                              |----- RT    => nonRT # start_port_dist/end_port_dist too large
        #$s_type = session_type() ;
        #chomp $s_type ;
        } elsif ($main::p_feed_pars_num{$id} == 2 && $main::p_sig_name{$id} ne "NA" && ($s_type eq 'flat' || $s_type eq 'ipo')) {
            $main::p_max_feed_pars{$main::p_sig_name{$id}}        = $main::p_feed_pars{$id} ;
            $main::p_max_man_dist{$main::p_sig_name{$id}}         = $main::p_man_dist{$id} ;
            $main::p_max_ideal_dist{$main::p_sig_name{$id}}       = $main::p_ideal_dist{$id} ;
            $main::p_max_real_dist{$main::p_sig_name{$id}}        = $main::p_real_dist{$id} ;
            $main::p_max_end_clk{$main::p_sig_name{$id}}          = $capt_clk ; 
            $main::p_start_routeRule{$id}                         = "NA" ; 
            $main::p_end_routeRule{$id}                           = "NA" ; 
            $main::p_sig_s_routeRule{$main::p_sig_name{$id}}      = $main::p_start_routeRule{$id} ; 
            $main::p_sig_e_routeRule{$main::p_sig_name{$id}}      = $main::p_end_routeRule{$id} ; 

            if ($start_unit !~ /_retime_partition_/ && $end_unit !~ /_retime_partition_/) {
                if ($main::p_is_mcp{$id}) {
                    $main::p_AdjParNRT2NRT_mcp{$main::p_sig_name{$id}} = $curr_vio ;
                } else {
                    $main::p_AdjParNRT2NRT_nonmcp{$main::p_sig_name{$id}} = $curr_vio ;
                }
            } elsif ($start_unit =~ /_retime_partition_/ && $end_unit !~ /_retime_partition_/) {
                $main::p_start_routeRule{$id}                    = GetStartRouteRule($curr_vio) ;
                $main::p_sig_s_routeRule{$main::p_sig_name{$id}} = $main::p_start_routeRule{$id} ; 
                $main::p_sig_is_mcp{$main::p_sig_name{$id}}      = $main::p_is_mcp{$id} ;
                if ($main::p_is_mcp{$id}) {
                    $main::p_AdjParRT2NRT_mcp{$main::p_sig_name{$id}} = $curr_vio ; 
                } else {
                    $main::p_AdjParRT2NRT_nonmcp{$main::p_sig_name{$id}} = $curr_vio ;
                }
            } elsif ($start_unit !~ /_retime_partition_/ && $end_unit =~ /_retime_partition_/) {
                $main::p_end_routeRule{$id}                      = GetEndRouteRule($curr_vio) ;
                $main::p_sig_e_routeRule{$main::p_sig_name{$id}} = $main::p_end_routeRule{$id} ; 
                $main::p_sig_is_mcp{$main::p_sig_name{$id}}      = $main::p_is_mcp{$id} ;
                if ($main::p_is_mcp{$id}) {
                    $main::p_AdjParNRT2RT_mcp{$main::p_sig_name{$id}} = $curr_vio ;
                } else {
                    $main::p_AdjParNRT2RT_nonmcp{$main::p_sig_name{$id}} = $curr_vio ;
                }
            } else {
                $main::p_start_routeRule{$id}                    = GetStartRouteRule($curr_vio) ;
                $main::p_end_routeRule{$id}                      = GetEndRouteRule($curr_vio) ;
                $main::p_sig_s_routeRule{$main::p_sig_name{$id}} = $main::p_start_routeRule{$id} ; 
                $main::p_sig_e_routeRule{$main::p_sig_name{$id}} = $main::p_end_routeRule{$id} ; 
                $main::p_sig_is_mcp{$main::p_sig_name{$id}}      = $main::p_is_mcp{$id} ;
                if ($main::p_is_mcp{$id}) {
                    $main::p_AdjParRT2RT_mcp{$main::p_sig_name{$id}} = $curr_vio ;
                } else {
                    $main::p_AdjParRT2RT_nonmcp{$main::p_sig_name{$id}} = $curr_vio ; 
                }   
            } 
        ##### intra-partition paths 
        } elsif ($main::p_feed_pars_num{$id} == 1 && $is_io == 0) {
            $main::p_intra_par_vio{$id}   = $main::p_ideal_dist{$id} ;
            $main::p_start_routeRule{$id} = GetStartRouteRule($curr_vio) ;
            $main::p_end_routeRule{$id}   = GetEndRouteRule($curr_vio) ;
        } elsif ($main::p_feed_pars_num{$id} == 1 && $is_io == 1) {
            $main::p_io_vio{$id} = $main::p_man_dist{$id} ;
            $main::p_start_routeRule{$id} = GetStartRouteRule($curr_vio) ;
            $main::p_end_routeRule{$id}   = GetEndRouteRule($curr_vio) ;
            if (exists $main::p_io_max_vio{$main::p_sig_name{$id}}) {
                if ($main::p_io_vio{$id} > $main::p_io_max_vio{$main::p_sig_name{$id}}) {
                    $main::p_io_max_vio{$main::p_sig_name{$id}} = $main::p_io_vio{$id} ;
                    $main::p_io_worst_id{$main::p_sig_name{$id}} = $id ;
                } 
            } else {
                    $main::p_io_max_vio{$main::p_sig_name{$id}}  = $main::p_io_vio{$id} ;
                    $main::p_io_worst_id{$main::p_sig_name{$id}} = $id ;
                }
        } else {
            if ($debug) {
                print $OUT "Not Rpted : $id\n" ;
            }
        }
    }
}

sub dump_reports {
    my ($top,$rep_name) = @_;

    #$s_type = session_type() ;
    #chomp $s_type ;

    open DistRep, ">${rep_name}.distance" or die $!;
    
    print "\n    Start to Dump Dist Rep              @ " . `date` ;

    print DistRep "#\n";
    print DistRep "# Manhattan: get_distance of ONLY startpoint, endpoint\n";
    print DistRep "# Ideal:     get_distance of ONLY startpoint, partition pins, endpoint\n";
    print DistRep "# Real:      get_distance of ALL (startpoint, combinational cells, partition pins, endpoint)\n";
    print DistRep "# Detour:    div(ideal / manhattan)..use this to metric identify bad partition pins\n";
    print DistRep "#\n";

    foreach my $id (sort keys %main::p_dist_rep) {
        print DistRep "$main::p_dist_rep{$id}\n" ;
    }

    close DistRep ;
    
    open IODistRep, ">${rep_name}.io.distance" or die $!;
    print "    Start to Dump IODist Rep            @ " . `date` ; 

    print IODistRep "#\n";
    print IODistRep "# Manhattan: get_distance of ONLY startpoint, endpoint\n";
    print IODistRep "# Ideal:     get_distance of ONLY startpoint, partition pins, endpoint\n";
    print IODistRep "# Real:      get_distance of ALL (startpoint, combinational cells, partition pins, endpoint)\n";
    print IODistRep "# Detour:    div(ideal / manhattan)..use this to metric identify bad partition pins\n";
    print IODistRep "#\n";

    foreach my $id (sort keys %main::p_io_dist_rep) {
        print IODistRep "$main::p_io_dist_rep{$id}\n" ;
    }
    
    close IODistRep ;
    
    open MannPlot, ">${rep_name}.manhattan_plot.medic" or die $!;
    print "    Start to Dump MannPlot Rep          @ " . `date` ; 

    foreach my $id (sort keys %main::p_mann_plot_rep) {
        print MannPlot "$main::p_mann_plot_rep{$id}\n" ;
    }

    close MannPlot ;

    open RetimeReport,  ">${rep_name}.full_retime.report" or die $!;
    print "    Start to Dump Full Retime Rep       @ " . `date` ;

    print RetimeReport "# Each column is: longest-net-name, end-clk, feed-pars, mann-dist, ideal-dsit\n";
    print RetimeReport "# !!! ideal_dist and feed_pars could be inaccurate without <par>_hfp.def\n";
    print RetimeReport "# !!! Please pay attention to the path if it's been commented out \n\n";

    foreach my $id (sort keys %main::p_retime_rep) {
        print RetimeReport "$main::p_retime_rep{$id}\n" ;
    }    

    close RetimeReport ;

    open DetourPlot, ">${rep_name}.detour_plot.medic" or die $!;
    print "    Start to Dump Detour plot file      @ " . `date` ;
        
    foreach my $id (sort keys %main::p_detour_plot_rep) {
        print DetourPlot "$main::p_detour_plot_rep{$id}\n" ;
    }

    close DetourPlot ;

    #$s_type = session_type() ;
    #chomp $s_type ;

    open FeedthrParRep , "> ${rep_name}.uniq_retime_feedpar.report" or die $! ;;
    print "    Start to Dump Feedthr Pars file     @ " . `date` ;
    
    print FeedthrParRep "# The signals through partitions :\n\n" ;
    print FeedthrParRep "\n##################################\n" ;
    print FeedthrParRep "# nonRT => nonRT && IS_MCP == 0  #\n" ;
    print FeedthrParRep "##################################\n" ;

    foreach my $sig_name (sort keys %main::p_max_nonrt_nonmcp_worst_id) {
        my $capt_clk = $main::p_max_nonrt_nonmcp_end_clk{$sig_name} ;
        $main::p_max_nonrt_nonmcp_clks{$capt_clk} = $main::p_clk_period{$main::p_max_nonrt_nonmcp_worst_id{$sig_name}} ; 
        if ($s_type eq 'noscan' || $s_type eq 'feflat') {
            $main::p_FeedParNRT2NRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_nonmcp_man_dist{$sig_name} ;
        } elsif ($s_type eq 'flat') {
            $main::p_FeedParNRT2NRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_nonmcp_ideal_dist{$sig_name} ; 
        } elsif ($s_type eq 'ipo') {
            $main::p_FeedParNRT2NRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_nonmcp_real_dist{$sig_name} ;
        }
    }
    
    foreach my $end_clk (sort {$main::p_max_nonrt_nonmcp_clks{$a} <=> $main::p_max_nonrt_nonmcp_clks{$b}} keys %main::p_max_nonrt_nonmcp_clks) {
        print FeedthrParRep "\n# Clk : $end_clk  Period : $main::p_max_nonrt_nonmcp_clks{$end_clk}ns\n" ;
        printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "is_mcp: <\${is_mcp}>" ) ;  
        foreach my $sig_name (sort {$main::p_FeedParNRT2NRT_nonmcp{$end_clk}{$b} <=> $main::p_FeedParNRT2NRT_nonmcp{$end_clk}{$a}} keys %{$main::p_FeedParNRT2NRT_nonmcp{$end_clk}}) {
            printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%15s # nonRT  & MCP == 0\n", $sig_name, $main::p_max_nonrt_nonmcp_feed_pars{$sig_name}, $main::p_max_nonrt_nonmcp_end_clk{$sig_name}, $main::p_max_nonrt_nonmcp_man_dist{$sig_name}, $main::p_max_nonrt_nonmcp_ideal_dist{$sig_name}, $main::p_max_nonrt_nonmcp_real_dist{$sig_name}, $main::p_max_nonrt_nonmcp_worst_id{$sig_name}, "is_mcp: $main::p_nonrt_nonmcp_sig_is_mcp{$sig_name}") ;
        }
    }

    print FeedthrParRep "\n##################################\n" ;
    print FeedthrParRep "# nonRT => nonRT && IS_MCP == 1  #\n" ;
    print FeedthrParRep "##################################\n" ;

    foreach my $sig_name (sort keys %main::p_max_nonrt_mcp_worst_id) {
        my $capt_clk = $main::p_max_nonrt_mcp_end_clk{$sig_name} ;
        $main::p_max_nonrt_mcp_clks{$capt_clk} = $main::p_clk_period{$main::p_max_nonrt_mcp_worst_id{$sig_name}} ; 
        if ($s_type eq 'noscan' || $s_type eq 'feflat') {
            $main::p_FeedParNRT2NRT_mcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_mcp_man_dist{$sig_name} ;
        } elsif ($s_type eq 'flat') { 
            $main::p_FeedParNRT2NRT_mcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_mcp_ideal_dist{$sig_name} ;
        } elsif ($s_type eq 'ipo') {
            $main::p_FeedParNRT2NRT_mcp{$capt_clk}{$sig_name} = $main::p_max_nonrt_mcp_real_dist{$sig_name} ;
        }
    }

    foreach my $end_clk (sort {$main::p_max_nonrt_mcp_clks{$a} <=> $main::p_max_nonrt_mcp_clks{$b}} keys %main::p_max_nonrt_mcp_clks) {
        print FeedthrParRep "\n# Clk : $end_clk  Period : $main::p_max_nonrt_mcp_clks{$end_clk}\n" ;
        printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "is_mcp: <\${is_mcp}>" ) ;

        foreach my $sig_name (sort {$main::p_FeedParNRT2NRT_mcp{$end_clk}{$b} <=> $main::p_FeedParNRT2NRT_mcp{$end_clk}{$a}} keys %{$main::p_FeedParNRT2NRT_mcp{$end_clk}}) {
            printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%15s # nonRT  & MCP == 1\n", $sig_name, $main::p_max_nonrt_mcp_feed_pars{$sig_name}, $main::p_max_nonrt_mcp_end_clk{$sig_name}, $main::p_max_nonrt_mcp_man_dist{$sig_name}, $main::p_max_nonrt_mcp_ideal_dist{$sig_name}, $main::p_max_nonrt_mcp_real_dist{$sig_name}, $main::p_max_nonrt_mcp_worst_id{$sig_name}, "is_mcp: $main::p_nonrt_mcp_sig_is_mcp{$sig_name}") ;
        }
    }

    print FeedthrParRep "\n###########################################################\n" ;
    print FeedthrParRep "# (nonRT => RT || RT => nonRT || RT => RT) && IS_MCP == 0 #\n" ;
    print FeedthrParRep "###########################################################\n" ;
    

    foreach my $sig_name (sort keys %main::p_max_rt_nonmcp_worst_id) {
        my $capt_clk = $main::p_max_rt_nonmcp_end_clk{$sig_name} ;
        $main::p_max_rt_nonmcp_clks{$capt_clk} = $main::p_clk_period{$main::p_max_rt_nonmcp_worst_id{$sig_name}} ;
        if ($s_type eq 'noscan' || $s_type eq 'feflat') {
            $main::p_FeedParRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_rt_nonmcp_man_dist{$sig_name} ;
        } elsif ($s_type eq 'flat') {
            $main::p_FeedParRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_rt_nonmcp_ideal_dist{$sig_name} ;
        } elsif ($s_type eq 'ipo') {
            $main::p_FeedParRT_nonmcp{$capt_clk}{$sig_name} = $main::p_max_rt_nonmcp_real_dist{$sig_name} ;
        }
    }


    foreach my $end_clk (sort {$main::p_max_rt_nonmcp_clks{$a} <=> $main::p_max_rt_nonmcp_clks{$b}} keys %main::p_max_rt_nonmcp_clks) {
        print FeedthrParRep "\n# Clk : $end_clk  Period : $main::p_max_rt_nonmcp_clks{$end_clk}\n" ;
        printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;  

        foreach my $sig_name (sort {$main::p_FeedParRT_nonmcp{$end_clk}{$b} <=> $main::p_FeedParRT_nonmcp{$end_clk}{$a}} keys %{$main::p_FeedParRT_nonmcp{$end_clk}}) {
            printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s # RT  & MCP == 0\n", $sig_name, $main::p_max_rt_nonmcp_feed_pars{$sig_name}, $main::p_max_rt_nonmcp_end_clk{$sig_name}, $main::p_max_rt_nonmcp_man_dist{$sig_name}, $main::p_max_rt_nonmcp_ideal_dist{$sig_name}, $main::p_max_rt_nonmcp_real_dist{$sig_name}, $main::p_max_rt_nonmcp_worst_id{$sig_name}, $main::p_max_rt_nonmcp_s_routeRule{$sig_name}, $main::p_max_rt_nonmcp_e_routeRule{$sig_name}, "is_mcp: $main::p_rt_nonmcp_sig_is_mcp{$sig_name}") ;
        }
    }

    print FeedthrParRep "\n###########################################################\n" ;
    print FeedthrParRep "# (nonRT => RT || RT => nonRT || RT => RT) && IS_MCP == 1 #\n" ;
    print FeedthrParRep "###########################################################\n" ;

    foreach my $sig_name (sort keys %main::p_max_rt_mcp_worst_id) {
        my $capt_clk = $main::p_max_rt_mcp_end_clk{$sig_name} ;
        $main::p_max_rt_mcp_clks{$capt_clk} = $main::p_clk_period{$main::p_max_rt_mcp_worst_id{$sig_name}} ;
        if ($s_type eq 'noscan' || $s_type eq 'feflat') {
            $main::p_FeedParRT_mcp{$capt_clk}{$sig_name} = $main::p_max_rt_mcp_man_dist{$sig_name} ;
        } elsif ($s_type eq 'flat') {
            $main::p_FeedParRT_mcp{$capt_clk}{$sig_name} = $main::p_max_rt_mcp_ideal_dist{$sig_name} ;
        } elsif ($s_type eq 'ipo') {
            $main::p_FeedParRT_mcp{$capt_clk}{$sig_name} = $main::p_max_rt_mcp_real_dist{$sig_name} ;
        }
    }


    foreach my $end_clk (sort {$main::p_max_rt_mcp_clks{$a} <=> $main::p_max_rt_mcp_clks{$b}} keys %main::p_max_rt_mcp_clks) {
        print FeedthrParRep "\n# Clk : $end_clk  Period : $main::p_max_rt_mcp_clks{$end_clk}\n" ;
        printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;  

        foreach my $sig_name (sort {$main::p_FeedParRT_mcp{$end_clk}{$b} <=> $main::p_FeedParRT_mcp{$end_clk}{$a}} keys %{$main::p_FeedParRT_mcp{$end_clk}}) {
            printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s # RT  & MCP == 1\n", $sig_name, $main::p_max_rt_mcp_feed_pars{$sig_name}, $main::p_max_rt_mcp_end_clk{$sig_name}, $main::p_max_rt_mcp_man_dist{$sig_name}, $main::p_max_rt_mcp_ideal_dist{$sig_name}, $main::p_max_rt_mcp_real_dist{$sig_name}, $main::p_max_rt_mcp_worst_id{$sig_name}, $main::p_max_rt_mcp_s_routeRule{$sig_name}, $main::p_max_rt_mcp_e_routeRule{$sig_name}, "is_mcp: $main::p_rt_mcp_sig_is_mcp{$sig_name}") ;
        }
    }

    print FeedthrParRep "\n###############################\n" ;
    print FeedthrParRep "# Detour caused by comb logic #\n" ;
    print FeedthrParRep "###############################\n" ;

    printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

    foreach my $sig_name (sort keys %main::p_comb_detour_start_routeRule) {
        printf FeedthrParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", $sig_name, $main::p_comb_detour_feed_pars{$sig_name}, $main::p_comb_detour_end_clk{$sig_name}, $main::p_comb_detour_man_dist{$sig_name}, $main::p_comb_detour_ideal_dist{$sig_name}, $main::p_comb_detour_real_dist{$sig_name}, $main::p_comb_detour_worst_id{$sig_name}, $main::p_comb_detour_start_routeRule{$sig_name}, $main::p_comb_detour_end_routeRule{$sig_name}, "is_mcp: $main::p_comb_detour_is_mcp{$sig_name}") ;
    }

    close FeedthrParRep ;


    #$s_type = session_type() ;
    #chomp $s_type ;
    if ($s_type eq 'flat' || $s_type eq 'ipo') {
        print "    Start to Dump Neighbor Pars file    @ " . `date` ;
        open AdjacentParRep, "> ${rep_name}.uniq_retime_adjpar.report" or die $! ;;

        print AdjacentParRep "# This is for neighbor partitions reports only on flat/ipo NL. \n\n" ;

        print AdjacentParRep "\n#################################\n" ; 
        print AdjacentParRep "# nonRT => nonRT && IS_MCP == 0 #\n" ;
        print AdjacentParRep "#################################\n" ;

        my %adj_par_nrt2nrt_nonmcp_clk_period = () ;
        my %adj_par_nrt2nrt_nonmcp_dist       = () ;
        my %adj_par_nrt2nrt_nonmcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParNRT2NRT_nonmcp) {
            my $vio      = $main::p_AdjParNRT2NRT_nonmcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ; 
            $adj_par_nrt2nrt_nonmcp_clk_period{$end_clk} = $main::p_clk_period{$id} ; 
            if ($s_type eq 'ipo') {
                $adj_par_nrt2nrt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;        
            } else {
                $adj_par_nrt2nrt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name}    = GetSourcePortDist ($vio) ; 
            $main::p_max_dest_port_dist{$sig_name}      = GetDestPortDist ($vio) ;
            $adj_par_nrt2nrt_nonmcp_worst_id{$sig_name} = $id ;
        }

        foreach my $end_clk (sort {$adj_par_nrt2nrt_nonmcp_clk_period{$a} <=> $adj_par_nrt2nrt_nonmcp_clk_period{$b}} keys %adj_par_nrt2nrt_nonmcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_nrt2nrt_nonmcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_nrt2nrt_nonmcp_dist{$end_clk}{$b} <=> $adj_par_nrt2nrt_nonmcp_dist{$end_clk}{$a}} keys %{$adj_par_nrt2nrt_nonmcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%15s # nonRT => nonRT & MCP == 0\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_nrt2nrt_nonmcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, "is_mcp: $main::p_is_mcp{$adj_par_nrt2nrt_nonmcp_worst_id{$sig_name}}") ;
            }
        }

        print AdjacentParRep "\n#################################\n" ;
        print AdjacentParRep "# nonRT => nonRT && IS_MCP == 1 #\n" ;
        print AdjacentParRep "#################################\n" ;

        my %adj_par_nrt2nrt_mcp_clk_period = () ;
        my %adj_par_nrt2nrt_mcp_dist       = () ;
        my %adj_par_nrt2nrt_mcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParNRT2NRT_mcp) {
            my $vio      = $main::p_AdjParNRT2NRT_mcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_nrt2nrt_mcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_nrt2nrt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_nrt2nrt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name} = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}   = GetDestPortDist ($vio) ;
            $adj_par_nrt2nrt_mcp_worst_id{$sig_name} = $id ;
        }

        foreach my $end_clk (sort {$adj_par_nrt2nrt_mcp_clk_period{$a} <=> $adj_par_nrt2nrt_mcp_clk_period{$b}} keys %adj_par_nrt2nrt_mcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_nrt2nrt_mcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "is_mcp: <\${is_mcp}>" ) ;
            
            foreach my $sig_name (sort {$adj_par_nrt2nrt_mcp_dist{$end_clk}{$b} <=> $adj_par_nrt2nrt_mcp_dist{$end_clk}{$a}} keys %{$adj_par_nrt2nrt_mcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%15s # nonRT => nonRT & MCP == 1\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_nrt2nrt_mcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, "is_mcp: $main::p_is_mcp{$adj_par_nrt2nrt_mcp_worst_id{$sig_name}}") ;
            }
        }

        print AdjacentParRep "\n##############################\n" ;
        print AdjacentParRep "# RT => nonRT && IS_MCP == 0 #\n" ;
        print AdjacentParRep "##############################\n" ;

        my %adj_par_rt2nrt_nonmcp_clk_period = () ;
        my %adj_par_rt2nrt_nonmcp_dist       = () ;
        my %adj_par_rt2nrt_nonmcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParRT2NRT_nonmcp) {
            my $vio      = $main::p_AdjParRT2NRT_nonmcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_rt2nrt_nonmcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_rt2nrt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_rt2nrt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name}   = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}     = GetDestPortDist ($vio) ;
            $adj_par_rt2nrt_nonmcp_worst_id{$sig_name} = $id ;
        }

        foreach my $end_clk (sort {$adj_par_rt2nrt_nonmcp_clk_period{$a} <=> $adj_par_rt2nrt_nonmcp_clk_period{$b}} keys %adj_par_rt2nrt_nonmcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_rt2nrt_nonmcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_rt2nrt_nonmcp_dist{$end_clk}{$b} <=> $adj_par_rt2nrt_nonmcp_dist{$end_clk}{$a}} keys %{$adj_par_rt2nrt_nonmcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # RT => nonRT & MCP == 0\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_rt2nrt_nonmcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n##############################\n" ;
        print AdjacentParRep "# RT => nonRT && IS_MCP == 1 #\n" ;
        print AdjacentParRep "##############################\n" ;

        my %adj_par_rt2nrt_mcp_clk_period = () ;
        my %adj_par_rt2nrt_mcp_dist       = () ;
        my %adj_par_rt2nrt_mcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParRT2NRT_mcp) {
            my $vio      = $main::p_AdjParRT2NRT_mcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_rt2nrt_mcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_rt2nrt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_rt2nrt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name} = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}   = GetDestPortDist ($vio) ;
            $adj_par_rt2nrt_mcp_worst_id{$sig_name}  = $id ;
        }

        foreach my $end_clk (sort {$adj_par_rt2nrt_mcp_clk_period{$a} <=> $adj_par_rt2nrt_mcp_clk_period{$b}} keys %adj_par_rt2nrt_mcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_rt2nrt_mcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_rt2nrt_mcp_dist{$end_clk}{$b} <=> $adj_par_rt2nrt_mcp_dist{$end_clk}{$a}} keys %{$adj_par_rt2nrt_mcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # RT => nonRT & MCP == 1\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_rt2nrt_mcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n##############################\n" ;
        print AdjacentParRep "# nonRT => RT && IS_MCP == 0 #\n" ;
        print AdjacentParRep "##############################\n" ;

        my %adj_par_nrt2rt_nonmcp_clk_period = () ;
        my %adj_par_nrt2rt_nonmcp_dist       = () ;
        my %adj_par_nrt2rt_nonmcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParNRT2RT_nonmcp) {
            my $vio      = $main::p_AdjParNRT2RT_nonmcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_nrt2rt_nonmcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_nrt2rt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_nrt2rt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name}   = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}     = GetDestPortDist ($vio) ;
            $adj_par_nrt2rt_nonmcp_worst_id{$sig_name} = $id ;
        }

        foreach my $end_clk (sort {$adj_par_nrt2rt_nonmcp_clk_period{$a} <=> $adj_par_nrt2rt_nonmcp_clk_period{$b}} keys %adj_par_nrt2rt_nonmcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_nrt2rt_nonmcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_nrt2rt_nonmcp_dist{$end_clk}{$b} <=> $adj_par_nrt2rt_nonmcp_dist{$end_clk}{$a}} keys %{$adj_par_nrt2rt_nonmcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # nonRT => RT & MCP == 0\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_nrt2rt_nonmcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n##############################\n" ;
        print AdjacentParRep "# nonRT => RT && IS_MCP == 1 #\n" ;
        print AdjacentParRep "##############################\n" ;
    
        my %adj_par_nrt2rt_mcp_clk_period = () ;
        my %adj_par_nrt2rt_mcp_dist       = () ;
        my %adj_par_nrt2rt_mcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParNRT2RT_mcp) {
            my $vio      = $main::p_AdjParNRT2RT_mcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_nrt2rt_mcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_nrt2rt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_nrt2rt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name} = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}   = GetDestPortDist ($vio) ;
            $adj_par_nrt2rt_mcp_worst_id{$sig_name}  = $id ;
        }

        foreach my $end_clk (sort {$adj_par_nrt2rt_mcp_clk_period{$a} <=> $adj_par_nrt2rt_mcp_clk_period{$b}} keys %adj_par_nrt2rt_mcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_nrt2rt_mcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;
            
            foreach my $sig_name (sort {$adj_par_nrt2rt_mcp_dist{$end_clk}{$b} <=> $adj_par_nrt2rt_mcp_dist{$end_clk}{$a}} keys %{$adj_par_nrt2rt_mcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # nonRT => RT & MCP == 1\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_nrt2rt_mcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n###########################\n" ;
        print AdjacentParRep "# RT => RT && IS_MCP == 0 #\n" ;
        print AdjacentParRep "###########################\n" ;

        my %adj_par_rt2rt_nonmcp_clk_period = () ;
        my %adj_par_rt2rt_nonmcp_dist       = () ;
        my %adj_par_rt2rt_nonmcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParRT2RT_nonmcp) {
            my $vio      = $main::p_AdjParRT2RT_nonmcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_rt2rt_nonmcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_rt2rt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_rt2rt_nonmcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name} = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}   = GetDestPortDist ($vio) ;
            $adj_par_rt2rt_nonmcp_worst_id{$sig_name}  = $id ;
        }

        foreach my $end_clk (sort {$adj_par_rt2rt_nonmcp_clk_period{$a} <=> $adj_par_rt2rt_nonmcp_clk_period{$b}} keys %adj_par_rt2rt_nonmcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_rt2rt_nonmcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_rt2rt_nonmcp_dist{$end_clk}{$b} <=> $adj_par_rt2rt_nonmcp_dist{$end_clk}{$a}} keys %{$adj_par_rt2rt_nonmcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # RT => RT & MCP == 0\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_real_dist{$sig_name}, $adj_par_rt2rt_nonmcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n###########################\n" ;
        print AdjacentParRep "# RT => RT && IS_MCP == 1 #\n" ;
        print AdjacentParRep "###########################\n" ;

        my %adj_par_rt2rt_mcp_clk_period = () ;
        my %adj_par_rt2rt_mcp_dist       = () ;
        my %adj_par_rt2rt_mcp_worst_id   = () ;
        foreach my $sig_name (sort keys %main::p_AdjParRT2RT_mcp) {
            my $vio      = $main::p_AdjParRT2RT_mcp{$sig_name} ;
            my $id       = attr_of_vio ('id' => $vio) ;
            my $end_clk  = $main::p_max_end_clk{$sig_name} ;
            $adj_par_rt2rt_mcp_clk_period{$end_clk} = $main::p_clk_period{$id} ;
            if ($s_type eq 'ipo') {
                $adj_par_rt2rt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_real_dist{$sig_name} ;
            } else {
                $adj_par_rt2rt_mcp_dist{$end_clk}{$sig_name} = $main::p_max_ideal_dist{$sig_name} ;
            }
            $main::p_max_source_port_dist{$sig_name} = GetSourcePortDist ($vio) ;
            $main::p_max_dest_port_dist{$sig_name}   = GetDestPortDist ($vio) ;
            $adj_par_rt2rt_mcp_worst_id{$sig_name}  = $id ;
        }

        foreach my $end_clk (sort {$adj_par_rt2rt_mcp_clk_period{$a} <=> $adj_par_rt2rt_mcp_clk_period{$b}} keys %adj_par_rt2rt_mcp_clk_period) {
            print AdjacentParRep "\n# Clk : $end_clk Period : $adj_par_rt2rt_mcp_clk_period{$end_clk}\n" ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s\n", "#<Signal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Worst ID>", "<Source Port Dist>", "<Dest Port Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>" ) ;

            foreach my $sig_name (sort {$adj_par_rt2rt_mcp_dist{$end_clk}{$b} <=> $adj_par_rt2rt_mcp_dist{$end_clk}{$a}} keys %{$adj_par_rt2rt_mcp_dist{$end_clk}}) {
                printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-15s\t%-25s\t%-25s\t%-30s\t%-30s\t%15s # RT => RT & MCP == 1\n", $sig_name, $main::p_max_feed_pars{$sig_name}, $main::p_max_end_clk{$sig_name}, $main::p_max_man_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $main::p_max_ideal_dist{$sig_name}, $adj_par_rt2rt_mcp_worst_id{$sig_name}, $main::p_max_source_port_dist{$sig_name}, $main::p_max_dest_port_dist{$sig_name}, $main::p_sig_s_routeRule{$sig_name}, $main::p_sig_e_routeRule{$sig_name}, "is_mcp: $main::p_sig_is_mcp{$sig_name}") ;
            }
        }

        print AdjacentParRep "\n############################\n" ;
        print AdjacentParRep "# Intra-Partition RT Paths #\n" ;
        print AdjacentParRep "############################\n" ;

        printf AdjacentParRep ("%-20s\t%-20s\t%-20s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", "#<ID>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>") ;
        
        foreach my $id (sort {$main::p_intra_par_vio{$b} <=> $main::p_intra_par_vio{$a}} keys %main::p_intra_par_vio) {
            printf AdjacentParRep ("%-20s\t%-20s\t%-20s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", $id, $main::p_feed_pars{$id}, $main::p_end_clk{$id}, $main::p_man_dist{$id}, $main::p_ideal_dist{$id}, $main::p_real_dist{$id}, $main::p_start_routeRule{$id}, $main::p_end_routeRule{$id}, "is_mcp: $main::p_is_mcp{$id}") ;
        } 

        print AdjacentParRep "\n############\n" ;
        print AdjacentParRep "# IO paths #\n" ;
        print AdjacentParRep "############\n" ;
    
        printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", "#<Siganal Name>", "<Feedthr Partitions>", "<End Clk>", "<Man Dist>", "<Ideal Dist>", "<Real Dist>", "<Start routeRule>", "<End routeRule>", "is_mcp: <\${is_mcp}>") ;
    
        foreach my $sig_name (sort {$main::p_io_max_vio{$b} <=> $main::p_io_max_vio{$a}} keys %main::p_io_max_vio) {
            my $id = $main::p_io_worst_id{$sig_name} ;
            printf AdjacentParRep ("%-80s\t%-40s\t\t%-20s\t%-15s\t%-15s\t%-15s\t%-30s\t%-30s\t%15s\n", $sig_name, $main::p_feed_pars{$id}, $main::p_end_clk{$id}, $main::p_man_dist{$id}, $main::p_ideal_dist{$id}, $main::p_real_dist{$id}, $main::p_start_routeRule{$id}, $main::p_end_routeRule{$id}, "is_mcp: $main::p_is_mcp{$id}") ;
        }
    
        close AdjacentParRep ;
    }

    print "    Start to Dump Histogram file        @ " . `date` ;
    foreach my $vio (@main::selected_vios) {
        my $id = attr_of_vio (id => $vio) ;
        #print "$id m: $main::p_man_dist{$id} r: $main::p_real_dist{$id} i: $main::p_ideal_dist{$id} s: $main::share_variables{"DIST_HISTOGRAM_SCALE"} \n" ;
        $main::p_dist_histogram_mann{ get_histogram_idx($main::p_man_dist{$id}, $main::share_variables{"DIST_HISTOGRAM_SCALE"}, $main::share_variables{"DIST_HISTOGRAM_MIN"}, $main::share_variables{"DIST_HISTOGRAM_MAX"}) } += 1;
        $main::p_dist_histogram_real{ get_histogram_idx($main::p_real_dist{$id}, $main::share_variables{"DIST_HISTOGRAM_SCALE"}, $main::share_variables{"DIST_HISTOGRAM_MIN"}, $main::share_variables{"DIST_HISTOGRAM_MAX"}) } += 1;
        $main::p_dist_histogram_ideal{ get_histogram_idx($main::p_ideal_dist{$id}, $main::share_variables{"DIST_HISTOGRAM_SCALE"}, $main::share_variables{"DIST_HISTOGRAM_MIN"}, $main::share_variables{"DIST_HISTOGRAM_MAX"}) } += 1;
        $main::p_detour_histogram{ get_histogram_idx_decimal($main::p_detour_ratio{$id}, 0.1, 1.0, 2.0) } += 1;
        # re-annotate the vio attribute
        set_vio_attr($vio, man_distance           => $main::p_man_dist{$id}); 
        set_vio_attr($vio, detour_ratio           => $main::p_detour_ratio{$id}) ;
        set_vio_attr($vio, feed_pars              => $main::p_feed_pars{$id});
        set_vio_attr($vio, feed_pars_num          => $main::p_feed_pars_num{$id});
        set_vio_attr($vio, par_num                => $main::p_par_num{$id});
        set_vio_attr($vio, real_distance          => $main::p_real_dist{$id});
        set_vio_attr($vio, ideal_distance         => $main::p_ideal_dist{$id});
        set_vio_attr($vio, is_mcp                 => $main::p_is_mcp{$id});
        set_vio_attr($vio, mcp_setup_num          => $main::p_mcp_setup_num{$id});
        set_vio_attr($vio, bin_man_distance       => $main::p_bin_man_distance{$id});
        set_vio_attr($vio, bin_ideal_distance     => $main::p_bin_ideal_distance{$id});
        set_vio_attr($vio, bin_real_distance      => $main::p_bin_real_distance{$id});
        set_vio_attr($vio, sig_name               => $main::p_sig_name{$id});
    }

    open Histogram,     ">${rep_name}.distance.HISTOGRAM" or die $!;

    printf Histogram ("%s\n\n", "Summary for ${top}, report=${rep_name}") ;
    printf Histogram ("%7s\t%-10s\t%-10s\t%-10s\n", " ", "Mann", "Ideal", "Real") ;
    ## syntax is Numeric sort
    foreach my $loop (sort {$a <=> $b} (keys %main::p_dist_histogram_mann)) {
        my $sign = " ";
        if ($loop eq $main::share_variables{"DIST_HISTOGRAM_MIN"}) { $sign = "<"; }
        if ($loop eq $main::share_variables{"DIST_HISTOGRAM_MAX"}) { $sign = ">"; }
        printf Histogram ("%s %5s\t%-10s\t%-10s\t%-10s\n", $sign, $loop, $main::p_dist_histogram_mann{$loop}, $main::p_dist_histogram_ideal{$loop}, $main::p_dist_histogram_real{$loop}) ;
    }

    print Histogram "\n\tMultiplier\n";
    ## syntax is Numeric sort
    foreach my $loop (sort {$a <=> $b} (keys %main::p_detour_histogram)) {
        my $sign = " ";
        if ($loop == 2.0) { $sign = ">"; }
        printf Histogram ("%s %5s\t%-10s\n", $sign, $loop, $main::p_detour_histogram{$loop}) ;
    }

    close (Histogram);

    print "    Start to Retime/Detour Summary file @ " . `date` . "\n" ;

    open RDSum,         ">${rep_name}.retime_detour.sum" or die $!;

    my $RETIME_DISTANCE                 = $main::share_variables{"RETIME_DISTANCE"};
    my $RETIME_DISTANCE_for_inter       = (int (($RETIME_DISTANCE/200) * 0.8)) * 100;
    my $DIST_HISTOGRAM_SCALE            = $main::share_variables{"DIST_HISTOGRAM_SCALE"};
    my $DIST_HISTOGRAM_MAX              = $main::share_variables{"DIST_HISTOGRAM_MAX"};
    my $MAN_DISTANCE_max_step           = (int (($DIST_HISTOGRAM_MAX - $RETIME_DISTANCE)/$DIST_HISTOGRAM_SCALE) + 1 );
    my $MAN_DISTANCE_max_step_for_inter = (int (($DIST_HISTOGRAM_MAX/2 - $RETIME_DISTANCE_for_inter)/$DIST_HISTOGRAM_SCALE) + 1);

    print RDSum "For this summary, per stage retime distance is $RETIME_DISTANCE for intra-chiplet, and it is $RETIME_DISTANCE_for_inter for inter-chiplet retime analysis.\n\n";

    print RDSum "Summary for ${top} inter-chiplet distance based on report ${rep_name}: \n";
    print RDSum (join "\n", (report_vios (-filter=>"is_io", -by=>"end_clk", -show=>"bin(man_distance,-step, ${RETIME_DISTANCE_for_inter}, ${DIST_HISTOGRAM_SCALE}, ${MAN_DISTANCE_max_step_for_inter}) worst(man_distance) count(id)"))) ;
    print RDSum "\n\n" ;
    
    print RDSum "Summary for ${top} intra-chiplet distance based on report ${rep_name}: \n";

    print RDSum (join "\n", (report_vios (-filter=>"!is_io", -by=>"end_clk", -show=>"bin(man_distance,-step, ${RETIME_DISTANCE}, ${DIST_HISTOGRAM_SCALE}, ${MAN_DISTANCE_max_step}) worst(man_distance) count(id)"))) ;
    print RDSum "\n\n" ;

    print RDSum "Summary for ${top} inter-chiplet detour based on report ${rep_name}: \n";

    print RDSum (join "\n", (report_vios (-filter=>"is_io and 'man_distance > ${RETIME_DISTANCE}' and 'detour_ratio >= 1.2'", -by=>"end_clk", -show=>"bin(detour_ratio,-step, 1.2,0.1,9) worst(detour_ratio) count(id)"))) ;
    print RDSum "\n\n" ;
    
    print RDSum "Summary for ${top} intra-chiplet detour based on report ${rep_name}: \n";

    print RDSum (join "\n", (report_vios (-filter=>"!is_io and 'man_distance > ${RETIME_DISTANCE}' and 'detour_ratio >= 1.2'", -by=>"end_clk", -show=>"bin(detour_ratio,-step, 1.2,0.1,9) worst(detour_ratio) count(id)"))) ; 

    close (RDSum);
    
    #print "    Start to unified retime report file @ " . `date` . "\n" ;
    #uniqReport($rep_name);
    
    print "    All the processes are done done     @ " . `date` ;
}


sub session_type {
    my $block=get_root_parent();
    my $gv_file = $M_file_gv{$block};
    my ($view) = $gv_file =~ /.*\/$block\.(.*)\.gv/;
    # view = (noscan|mbist|feflat|flat|layout/ipo/anno|pretp)
    if ($view =~ /noscan/ && $view !~ /pnr/) { #noscan.flat, noscan.par
            $view = "noscan";
    } elsif ($view =~ /ipo\d+$/) {
            $view = "anno";
    } elsif ($view =~ /FE_flat$/) {
            $view = "feflat";
    } elsif ($view =~ /flat$/) {
            $view = "flat";
    } elsif ($view =~ /pretp$/) {
            $view = "pretp";
    } elsif ($view =~ /mbist$/) {
            $view = "mbist";
    } elsif ($view =~ /noscan/ && $view =~ /pnr/) {
            $view = "noscan_pnrprecheck";
    }else {
            $view = "layout";
    }
    # there's issue if the top netlist is a soft link in TOT
    # example in GA103, nv_top.flat.gv.gz -> nv_top.FE_flat.gv.gz, the flow returns feflat instead of flat
    if ($ENV{TS_VIEW} ne "") {
        $type = $ENV{TS_VIEW};
    } else {
        $type = $view;
    }
    return ($type);
}

sub genUnitNameHash {
    my $config = CadConfig::factory();
    my %unit_hier_mapping;
    my @all_units = sort(keys(%{$config->{partitioning}->{units}}));
    my %hash_units;
    @all_units = grep { ++$hash_units{$_} < 2 } @all_units;

    ### Generate unit -> hier_name hash
    foreach my $unit (@all_units) {
        if (%{$config->{partitioning}->{units}{$unit}{"partition"}}) {
            foreach my $s_par (keys(%{$config->{partitioning}->{units}{$unit}{"partition"}})) {
                my @inst_names = split(',',$config->{partitioning}->{units}{$unit}{"partition"}{$s_par});
                foreach my $inst_name (@inst_names) {
                    $unit_hier_mapping{"$s_par/$inst_name"} = $unit;
                }
            }
        }
    }    
    return %unit_hier_mapping;
}

Tsub mapPinUnit => << 'END' ; 
    DESC {
        to get the unit for pin
        Usage : mapPinUnit <Pin Name>
    }
    ARGS {
        $pin
    }

    my @units = sort keys %main::unit_hier_mapping ;
 
    if ($#units == -1) {
        %main::unit_hier_mapping = genUnitNameHash ; 
    }

    foreach my $unit_hier (keys %main::unit_hier_mapping) {
        if ($pin =~ /$unit_hier/) {
            return $unit_hier_mapping{$unit_hier};
        }
    }
    return "";

END

Tsub load_port_map_file => << 'END' ;
    DESC {
        to load the port mapping file.
    }    
    ARGS {
    }    

    %M_noscan_port_mapping = () ; 

    my $proj = $ENV{NV_PROJECT} ;
    my $rev  = $ENV{USE_LAYOUT_REV} ;
    my %all_mods = map  ({$_ => 1} (get_modules "*")) ;
    my @chiplets = grep ((exists $all_mods{$_}), (all_chiplets)) ;

    foreach my $top (@chiplets) {
        my $ipo_dir      = $ENV{IPO_DIR} ;
        my $portmap_file = "${ipo_dir}/${top}/noscan_cfg/${top}.noscan.portmap" ;

        if (-e $portmap_file) {
            print "Loading $portmap_file ...\n" ;
            open IN, "$portmap_file" ;
            while (<IN>) {
                chomp ;
                my $line = $_ ; 
                if ($line =~ /^\w+/) {
                    # partition_module partition_port  cell cell_inst cell_pin
                    my ($par_ref, $par_port, $unit_name, $unit_inst, $unit_port) = split (" ", $line) ;
                    $par_port =~ s/^\\// ;
                    $M_noscan_port_mapping{$par_ref}{$par_port} = 1 ;
                }
            }
            close IN ; 
        } else {
            error "Can't find portmap file : $portmap_file\n" ;
        }
    }    

    return 1 ;

END

Tsub load_routeRules_files => << 'END' ;
    DESC  {
        to parse the route Rule files in mender session
    }
    ARGS {
    }

    %M_routeRules      = () ;
    %M_routeRules_pipe = () ;
    $chiplet_uc        = "" ;

    my %all_mods = map  ({$_ => 1} (get_modules "*")) ;
    my @chiplets = grep ((exists $all_mods{$_}), (all_chiplets)) ;

    my $litter = $CONFIG->{LITTER_NAME} ;
    my $chip_root = `depth` ;
    chomp $chip_root ;

    my $routeRulesdir = "$chip_root/ip/retime/retime/1.0/vmod/include/interface_retime" ;
    foreach my $block (@chiplets) {
        my $chiplet        = $block ;
        $chiplet           =~ s/NV_(.*)/$1/ ;
        $chiplet_uc     = uc $chiplet ;
        my $routeRulesFile = "$routeRulesdir/interface_retime_${litter}_${chiplet_uc}_routeRules.pm" ;
        if (-e $routeRulesFile) {
            `p4 sync $routeRulesFile` ;
            #print "Loading routeRules file : $routeRulesFile\n" ;
            load "$routeRulesFile" ;
        } else {
            next ;
        }
    }

    foreach my $chiplet (sort keys %M_routeRules) {
        foreach my $rule_name (sort keys %{$M_routeRules{$chiplet}}) {
            if (exists $M_routeRules{$chiplet}{$rule_name}{pipeline_steps}) {
                foreach my $pipe (split (",", $M_routeRules{$chiplet}{$rule_name}{pipeline_steps})) {
                    $M_routeRules_pipe{$chiplet}{$rule_name}{$pipe} = 1 ;
                }
            }
            if (exists $M_routeRules{$chiplet}{$rule_name}{tap}) {
                foreach my $tap_num (sort keys %{$M_routeRules{$chiplet}{$rule_name}{tap}}) {
                    if (exists $M_routeRules{$chiplet}{$rule_name}{tap}{$tap_num}{pipeline_steps}) {
                        foreach my $pipe (split (",", $M_routeRules{$chiplet}{$rule_name}{tap}{$tap_num}{pipeline_steps})) {
                            $M_routeRules_pipe{$chiplet}{$rule_name}{$pipe} = 1 ;
                        }
                    }
                }
            }
        }
    }

    return 1 ;

END

Tsub AddRouteRule => << 'END' ;
    DESC {
        a dummy function for parsing the interface_retime_*.pm files
    }
    ARGS {
        @args
    }

    my %in = @args ;

    my $rule_name = $in{name} ;
    foreach my $key (sort keys %in) {
        $M_routeRules{$chiplet_uc}{$rule_name}{$key} = $in{$key} ;
    }

    return 1 ;

END


Tsub dump_mcp_tcl_file => << 'END' ;
    DESC {
        not dump the tcl file by default. 
    }
    ARGS {
        $rep_name
    }

    open McpReport, ">${rep_name}.temp_mcp_for_missing_retime.tcl" or die $!;
    print "    Start to Dump MCP tcl file @ " . `date` . "\n" ;
    
    foreach my $id (sort keys %main::p_mcp_rep) {
        print McpReport "$main::p_mcp_rep{$id}\n" ; 
    }
    close McpReport ;
END

sub get_all_par_insts {
    my @all_par_insts = () ;
    my %allModules    = map  ({$_ => 1} (get_modules ("*"))) ;
    my @partRefs      = grep (exists $allModules{$_}, (all_partitions)) ; 
    foreach my $part (@partRefs) {
        push @all_par_insts, (get_cells_of $part) ;
    }
    return @all_par_insts ;
}

Tsub GetVioSigName => << 'END' ;
    DESC {
        to get the RTL net name of one violation.
    }
    ARGS {
        $vio
    }

    my $rtn = "NA" ;
    my @hier_pins = attr_of_vio ('inter_pin_array' => $vio) ;

    foreach my $pin (@hier_pins) {
        if (is_port $pin) {
            #push @sig_name, $pin ;
            return $pin ;
        } else {
            if (get_pin (-quiet => $pin)) {
                my ($par, $par_ref, $pin_name) = get_pin_context $pin ;
                if ((is_partition_module $par_ref) && (exists $M_noscan_port_mapping{$par_ref}{$pin_name})){
                    $pin_name =~ s/\[\d+\]$// ;
                    #push @sig_names, $pin_name ;
                    return $pin_name ;
                } else {
                    next ;
                }
            } else {
                next ;
            }
        }
    }

    return $rtn ;

END

Tsub GetSourcePortDist => << 'END' ;
    DESC {
        get the dist between violation start_pin and start_par_cell port.
    }
    ARGS {
        $vio
    }

    my $start_pin    = attr_of_vio (start_pin => $vio) ;
    my $start_par    = attr_of_vio (start_par => $vio) ;
    my @hier_port    = grep (attr_of_pin (is_hier, $_), split (" ", attr_of_vio (pin_list => $vio))) ;
    my @source_ports = () ;
    my $dist ;

    foreach my $port (@hier_port) {
        if (is_port $port) {
            push @source_ports, $port ;
        } else {
            my @pin_context = get_pin_context ($port) ;
            if ($pin_context[1] eq $start_par) {
                #$dist = get_dist ($port, $start_pin) ; 
                push @source_ports, $port ;
            }
        }
    }

    if ($#source_ports > 0) {
        if (is_port $source_ports[0]) {
            $dist = get_dist ($start_pin, $source_ports[1]) ;
        } else {
            $dist = get_dist ($start_pin, $source_ports[0]) ;
        }
    } elsif ($#source_ports == -1) {
        $dist = 0 ;
    } else {
        $dist = get_dist ($start_pin, $source_ports[0]) ;
    }

    return $dist ;

END

Tsub GetDestPortDist => << 'END' ;
    DESC {
        get the dist between violation end_pin and end_par_cell port.
    }
    ARGS {
        $vio
    }

    my $end_pin    = attr_of_vio (end_pin => $vio) ;
    my $end_par    = attr_of_vio (end_par => $vio) ;
    my @hier_port  = grep (attr_of_pin (is_hier, $_), split (" ", attr_of_vio (pin_list => $vio))) ;
    my @dest_ports = () ;
    my $dist ;

    foreach my $port (@hier_port) {
        if (is_port $port) {
            push @dest_ports, $port ;
        } else {
            my @pin_context = get_pin_context ($port) ;
            if ($pin_context[1] eq $end_par) {
                #$dist = get_dist ($port, $end_pin) ;
                push @dest_ports, $port ;
            }
        }
    }

    if ($#dest_ports > 0) {
        $dist = get_dist ($end_pin, $dest_ports[-1]) ;
    } elsif ($#dest_ports == -1) {
        $dist = 0 ;
    } else {
        $dist = get_dist ($end_pin, $dest_ports[0]) ;
    }

    return $dist ;

END

Tsub GetSourceCoor => << 'END' ;
    DESC {
        get the violation start pin coordinate.
    }
    ARGS {
        $vio
    }

    my $start_pin = attr_of_vio (start_pin => $vio) ;
    my @coor      = get_pin_xy $start_pin ;

    return "$coor[0],$coor[1]" ;

END

Tsub GetDestCoor => << 'END' ;
    DESC {
        get the violation end pin coordinate.
    }
    ARGS {
        $vio
    }

    my $end_pin = attr_of_vio (end_pin => $vio) ;
    my @coor    = get_pin_xy $end_pin ;

    return "$coor[0],$coor[1]" ;

END

Tsub GetStartRouteRule => << 'END' ;
    DESC {
        get the violation start pin routeRules.
    }
    ARGS {
        $vio
    }

    #my $start_pin = attr_of_vio (start_pin => $vio) ;
    my @pin_list  = split (" ", attr_of_vio (pin_list => $vio)) ;;
    my $start_pin = $pin_list[1] ;
    my $routeRule = get_rule_of_pin ($start_pin) ;

    return $routeRule ;

END

Tsub GetEndRouteRule => << 'END' ;
    DESC {
        get the violation end pin routeRules.
    }
    ARGS {
        $vio
    }

    my $end_pin = attr_of_vio (end_pin => $vio) ;
    my $routeRule = get_rule_of_pin ($end_pin) ;

    return $routeRule ;

END

Tsub get_rule_of_pin => << 'END';
    DESC {
        To get the routeRule of retime flop pin
    }
    ARGS {
        $pin   ,
    }

    if (!$pin) {
        error "No pin listed.\n" ;
        return () ;
    }

    if (get_port (-quiet, $pin)) {
        return "NA" ;
    }

    if (!(get_pins (-quiet, $pin))) {
        error "$pin not found.\n" ;
        return () ;
    }

    my $inst ;
    my $ref ;

    if (get_pin (-quiet, $pin)) {
        $inst = get_cell (-of => $pin) ;
        $ref  = get_ref ($inst) ;

        if (attr_of_ref ('is_merged_flop' => $ref)) {
            my $demerged_pin = get_demerged_name ($pin) ;
            $pin = $demerged_pin ;
            print "Demerged Pin : $pin\n" if (defined $d) ;
        }
    }

    if ($pin =~ /_retime_.*_RT/) {
            my $pipe = $pin ;
            $pipe =~ s/.*_RT.*?_(\S+?)\/.*\/.*/$1/ ;
            foreach my $chiplet (sort keys %M_routeRules_pipe) {
                foreach my $rule (sort keys %{$M_routeRules_pipe{$chiplet}}) {
                    if (exists $M_routeRules_pipe{$chiplet}{$rule}{$pipe}) {
                        #return "$chiplet $rule $pipe" ;
                        return "$rule" ;
                    }
                }
            }
            return "NA" ;
    } else {
        return "NA" ;
    }

END

Tsub load_def_region_files => << 'END' ;
    DESC {
        to load the def region files 
    }
    ARGS {
        -ipo_dir:$ipo_dir ,
        -legalize_region ,
    }

    set_cell_size_default (-macro => (10, 10)) ;
    set_rc_default_estimated ;
    set_xy_default (-centroid) ;

    if (!defined $opt_legalize_region) {
        set_eco_legal_placement (never) ;
    }

    my @loaded_files = get_files() ;
    my $rev          = $ENV{USE_LAYOUT_REV} ;
    my $proj         = $ENV{NV_PROJECT} ;
    my @files        = () ;

    
    foreach my $file (@loaded_files) {
        $file =~ s/\S+:\s+(\S+)\s*.*/$1/ ;
        push @files, $file ;
    }

    if (!defined $ipo_dir) {
        $ipo_dir = $ENV{IPO_DIR} ; 
    }

    my %all_mods   = map ({$_ => 1} (get_modules("*"))) ;
    my @macros     = grep (exists ($all_mods{$_}), (all_macros));
    my @parts      = grep (exists ($all_mods{$_}), (all_partitions));
    my @chiplets   = grep (exists ($all_mods{$_}), (all_chiplets));

    # load hcoff.data file ;
    my $hcoff_data = "${ipo_dir}/${proj}_top/control/hcoff.data";
    load "$hcoff_data" ;
    
    tcm_set_var ("TOOL(NAME)", dummy);

    my $type = session_type() ;

    
    ## split noscan/feflat with flat, flat load full def/hfp def. 
    if ( $type eq "noscan" | $type eq "feflat") {
        foreach my $part (@macros) {
    	    if (-e "${ipo_dir}/macros/${part}/netlists/${part}.def.gz") {
    	        load_once "${ipo_dir}/macros/${part}/netlists/${part}.def.gz" ;
    	    } elsif (-e "${ipo_dir}/macros/${part}/netlists/${part}.def") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.def" ;
    	    } elsif (-e "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def.gz") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def.gz" ;
    	    } elsif  (-e "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def" ;
    	    } else {
    	    	print "# no def file found for ${part}\n";
    	    }
        }
    
        foreach my $part (@parts, @chiplets) {
    	    if (-e "${ipo_dir}/${part}/control/${part}_fp.def") {
    	        load_once "${ipo_dir}/${part}/control/${part}_fp.def" ;
                my $retime_tcl = "${ipo_dir}/${part}/control/${part}_RETIME.tcl" ;
                my $cts_retime_tcl = "${ipo_dir}/${part}/control/${part}_nvq_region.cts_retime.tcl" ;
    	    	if (-e "$retime_tcl") {
                    if (grep ($_ eq $retime_tcl, @files)) {
                        print "Already loaded $retime_tcl\n" ;
                    } else {
                        load "$retime_tcl" ;
                    }
    	    	}
                if (-e "$cts_retime_tcl") {
                    if (grep ($_ eq $cts_retime_tcl, @files)) {
                        print "Already loaded $cts_retime_tcl\n" ;
                    } else {
                        load "$cts_retime_tcl" ;
                    }
                }

    	    	## custom timing region used by XBAR team
                my $cus_tcl = "${ipo_dir}/${part}/control/${part}_timing_region.tcl" ;
                if (-e "$cus_tcl") {
    	    		#if($type eq "noscan" && $part !~ /GAEX0/ && $part !~ /NV_gae_x0/) {
    	    		#	push_top($part); @list = ();
    	    		#	@list = get_cells_of("XTR_UNIT_WRA*");
    	    		#	push(@list, get_cells_of("DUMMY*"));
    	    		#	ungroup_inst(-flatten => @list);
    	    		#	pop_top();
    	    		#}
                    if (grep ($_ eq $cus_tcl, @files)) {
                        print "Already loaded $cus_tcl\n" ;   
                    } else {
                        load "$cus_tcl" ;
                    }
    	    	}
    	    } else {
    	    	print "# No fp_def found for ${part}\n"; 
    	    }
        }
    } elsif ( $type eq "flat" | $type eq "noscan_pnrprecheck" ) {
        foreach my $part (@macros) {
    	    if (-e "${ipo_dir}/macros/${part}/netlists/${part}.def.gz") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.def.gz" ;
    	    } elsif (-e "${ipo_dir}/macros/${part}/netlists/${part}.def") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.def" ;
    	    } elsif (-e "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def.gz") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def.gz" ;
    	    } elsif (-e "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def") {
    	    	load_once "${ipo_dir}/macros/${part}/netlists/${part}.FE_flat.def" ;
    	    } else {
    	    	print "# no def file found for ${part}\n";
    	    }
        }
    
        foreach my $part (@parts, @chiplets) {
            if (-e "${ipo_dir}/${part}/netlists/${part}.hfp.full.def.gz") {
            	load_once "${ipo_dir}/${part}/netlists/${part}.hfp.full.def.gz" ;
            } elsif (-e "${ipo_dir}/${part}/control/${part}_fp.def") {
            	load_once "${ipo_dir}/${part}/control/${part}_fp.def" ;
                my $retime_tcl = "${ipo_dir}/${part}/control/${part}_RETIME.tcl" ;
                my $cts_retime_tcl = "${ipo_dir}/${part}/control/${part}_nvq_region.cts_retime.tcl" ;
                if (-e "$retime_tcl") {
                    if (grep ($_ eq $retime_tcl, @files)) {
                        print "Already loaded $retime_tcl\n" ;
                    } else {
                        load "$retime_tcl" ;
                    }
                }
                if (-e "$cts_retime_tcl") {
                    if (grep ($_ eq $cts_retime_tcl, @files)) {
                        print "Already loaded $cts_retime_tcl\n" ;
                    } else {
                        load "$cts_retime_tcl" ;
                    }
                }

            	if (-e "${ipo_dir}/${part}/control/${part}.hfp.pins.def") {
            		load (-replace, "${ipo_dir}/${part}/control/${part}.hfp.pins.def") ;
            	}
    
            	## custom timing region used by XBAR team
                my $cus_tcl = "${ipo_dir}/${part}/control/${part}_timing_region.tcl" ;
                if (-e "$cus_tcl") {
                    if (grep ($_ eq $cus_tcl, @files)) {
                        print "Already loaded $cus_tcl\n" ;
                    } else {
                        load "$cus_tcl" ;
                    }
                } 
            } else {
            	print "# No def found for ${part}\n"; 
            }
        }
    }


    # to load the def medic for user/chiplet sepcial
    my $user = $ENV{USER} ;
    my $top  = get_top ;
    chomp $top ;

    if (-e "${proj}/timing_scripts/workflow/retime_detour_support/load_def.${user}.inc.medic") {
        load "${proj}/timing_scripts/workflow/retime_detour_support/load_def.${user}.inc.medic" ;
    }

    if (-e "${proj}/timing_scripts/workflow/retime_detour_support/load_def.${top}.inc.medic") {
        load "${proj}/timing_scripts/workflow/retime_detour_support/load_def.${top}.inc.medic" ;
    }

END

sub get_dist_pinArray {
   my @pin_list = @_;
   my $pin_list_size = scalar @pin_list;

   my $RC = 0;
   for ($idx = 0; $idx < $pin_list_size-1; $idx++) {
    $RC += get_dist (@pin_list[$idx] => @pin_list[$idx+1]);
   }

   return $RC;
}

sub get_histogram_idx {
   my ($distance , $scale, $min, $max) = @_;

   my $myIdx = int ($distance / $scale) * $scale ;
   if ( $myIdx > $max ) { $myIdx = $max; }
   if ( $myIdx < $min ) { $myIdx = $min; }

   return $myIdx;
}

sub get_histogram_idx_decimal {
   my ($distance , $scale, $min, $max) = @_;

   my $distance_mod = $distance * 10;
   my $scale_mod = $scale * 10;
   my $min_mod = $min * 10;
   my $max_mod = $max * 10;

   my $myIdx_mod = get_histogram_idx($distance_mod, $scale_mod, $min_mod, $max_mod);

   return ($myIdx_mod / 10);
}

Tsub get_vio_feeds => << 'END' ;
    DESC {
        to get the feedthr partitions for violation 
    }
    ARGS {
        $vio ,
    }
    
    my @rtn_pars = () ;

    my @vio_thr_pars = get_vio_thr_pars ($vio) ;
    
    if ($#vio_thr_pars == 0) {
        push @rtn_pars, $vio_thr_pars[0] ;
    } else {
        foreach my $i (1..$#vio_thr_pars) {
            my @thr_pars = @{$M_all_thr_pars{$vio_thr_pars[$i-1]}{$vio_thr_pars[$i]}} ;
            foreach my $thr_par (@thr_pars) {
                if ($thr_par ne $rtn_pars[-1]) {
                    push @rtn_pars, $thr_par ;
                } else {
                    next ;
                }
            }
        } 
    }
    
    return @rtn_pars ;

END

Tsub get_vio_thr_pars => << 'END' ;
    DESC {
        to get the violation partition insts 
    }
    ARGS {
        $vio
    }

    my @rtn = () ;

    my $top                = get_top() ;
    my $start_pin          = attr_of_vio ('start_pin' => $vio) ;
    my $end_pin            = attr_of_vio ('end_pin' => $vio) ;
    my $start_par_cell     = () ;
    my $end_par_cell       = () ;

    if (is_port $start_pin) {
        $start_par_cell = get_port_partition ($start_pin) ;    
    } else {
        $start_par_cell = attr_of_vio ('start_par_cell' => $vio) ;
    }    

    if (is_port $end_pin) {
        $end_par_cell = get_port_partition ($end_pin) ;
    } else {
        $end_par_cell = attr_of_vio ('end_par_cell' => $vio) ;
    }

    my @inter_pins = attr_of_vio ('inter_pin_array' => $vio) ;
    if ($#inter_pins == -1) {
        if ($start_par_cell eq $end_par_cell) {
            push @rtn, $start_par_cell ;
        } else {
            print "check path : $start_pin $end_pin\n" ;
        }
    } else {
        push @rtn, $start_par_cell ;
        foreach my $inter_pin (@inter_pins) {
            my ($par_inst, $par_ref, $par_pin) = get_pin_context ($inter_pin) ;
            if ((is_partition_module $par_ref) && ($par_inst ne $rtn[-1])) {
                push @rtn, $par_inst ;
            } elsif (is_chiplet_module $par_ref) {
                set_top $par_ref ;
                my $sub_par_inst = get_port_partition ($par_pin) ;
                $par_inst = "$par_inst/$sub_par_inst" ;
                set_top $top ;
                if ($par_inst ne $rtn[-1]) {
                    push @rtn, $par_inst ;
                }
            } else {
                next ;
            }
        }
        if ($end_par_cell ne $rtn[-1]) {
            push @rtn, $end_par_cell ;
        }
    }

    return @rtn ; 

END

Tsub get_all_par_thrs => << 'END' ;
    DESC {
        to get all the partition feedthr hash.
    }
    ARGS {
    }
        
    if ($#main::all_partitions =-1) {
        @main::all_partitions = get_all_par_insts () ;
    } 

    foreach my $par (@main::all_partitions) {
        get_neighbor_partitions ($par) ;
    }

    my $graph = Graph::Directed->new ;

    foreach my $par_inst (keys %M_all_neighbor_pars) {
        $graph->add_edge($par_inst, $_) for @{$M_all_neighbor_pars{$par_inst}};
    }

    my $APSP = $graph->APSP_Floyd_Warshall;

    foreach my $start_par ($APSP-> vertices ) {
        foreach my $end_par ($APSP-> vertices ) {
            my $pars_num = $APSP->path_length($start_par, $end_par);
            if ($pars_num) {
                my @thr_pars = $APSP->path_vertices($start_par, $end_par);
                @{$M_all_thr_pars{$start_par}{$end_par}} = @thr_pars ;
            } else {
                next ;
            }
        }
    }

END

Tsub get_neighbor_partitions => << 'END' ;
    DESC {
        internal function, to get all the neighor partitions ;
    }
    ARGS {
        $par_inst ,
    }

    my $top = get_top ;
    chomp $top ;
    set_chip_top $top ;
    set_tl_hier_abutment 1 ;

    my %par_bound     = () ;
    my @neighbor_pars = () ;
    my @all_n_pars    = tl_get_abutments $par_inst ;

    foreach my $par_info (@all_n_pars) {
        my $par_name   = $par_info->[0] ;
        my @sp_coor    = ($par_info->[2], $par_info->[3]) ;
        my @ep_coor    = ($par_info->[4], $par_info->[5]) ;
        my $bound_dist = get_dist (@sp_coor, @ep_coor) ;
        if (get_cells (-quiet => $par_name)) {
            if (exists $par_bound{$par_name}) {
                $par_bound{$par_name} = $par_bound{$par_name} + $bound_dist ;
            } else {
                $par_bound{$par_name} = $bound_dist ;
            }
        }
    }

    foreach my $par (sort keys %par_bound) {
        if ($par_bound{$par} > 200) {
            push @neighbor_pars, $par ;
        } else {
            # print "short bound : $par $par_inst $par_bound{$par}\n" ;
        }
    }

    @{$M_all_neighbor_pars{$par_inst}} = @neighbor_pars ;
    return ;

END

#sub get_feed_pars {
#    #give net, this routine gives list of pars that net has to cross
#    #make sure you run below commands before calling this function
#    #set_verilog_build_from_layout
#    #set_libs_default_none
#    #load_coff_data /home/gv100_layout/tot/layout/revP3.0/blocks/gv100_top/control/hcoff.data
#    #set_chip_top NV_gva_ff0
#    my $net = shift;
#    my $end_par_u = shift;
#    my $top = shift;
#    my $end_par = lc($end_par_u);
#    my @parList;
#    return if(is_power_net($net)); #skip power nets
#    my $driver; ($driver) = get_drivers ($net);
#    my @loads = get_loads($net);
#    return if(scalar(@loads) < 1);
#    if(is_port($driver)) {#if port assign one loads
#        ($driver) = shift(@loads);
#    }
#    my $dpar = partition_inst_of_pin ($driver);
#    my ($dx,$dy) = get_object_xy ($dpar);
#    foreach $load (@loads) { #loop through each load and get farthest load
#        next if(is_port($load)); #skip chiplet/nv_top ports
#        my @temp; #array to check farthest load
#        my $lpar = partition_inst_of_pin($load);
#        if($lpar ne $dpar) {
#           if ($end_par eq $lpar) {
#            if ( (($dpar !~ /TPC/) && ($lpar !~ /TPC/)) && (($dpar !~ /GAA0SC/) && ($lpar !~ /GAA0SC/))) {
#                if (!(tl_is_abutted ($lpar,$dpar))) {
#                my ($lx,$ly) = get_object_xy ($lpar);
#                my @flist = tl_get_crossed_edges(-min_cross => $dx,$dy,$lx,$ly,$dpar,$lpar);
#                  foreach $a (@flist) {
#                    my $par = $a->[0];#parInst  
#                    push(@temp, $par) if($temp[-1] ne $par);
#                  }
#                } else {
#                    push(@temp,$dpar);
#                    push(@temp,$lpar);
#                }
#            } else {
#                  push(@temp,$dpar);
#                  push(@temp,$lpar);
#            }
#          } else {
#              push(@temp,$dpar);
#              push(@temp,$end_par);
#          }
#    } else {
#        push(@temp,$lpar);
#    }
#   if(scalar(@temp) > scalar(@parList)) { #assing farthest load
#            @parList = @temp;
#        }
#  }
#    return @parList;
#}
#
#sub get_feed_pars_for_pars {
#    my $spar_u = shift;
#    my $spar = lc($spar_u);
#    my $epar_u = shift;
#    my $epar = lc($epar_u);
#    my @parList;
#    return $spar if($spar eq $epar);
#    if (!(tl_is_abutted ($spar,$epar))) {
#      my ($dx,$dy) = get_object_xy ($spar);
#      my ($lx,$ly) = get_object_xy ($epar);
#      my @flist = tl_get_crossed_edges(-min_cross => $dx,$dy,$lx,$ly,$spar,$epar);
#      foreach $a (@flist) {
#        my $par = $a->[0];#parInst  
#        push(@parList, $par) if($parList[-1] ne $par);
#      }
#    } else {
#          push(@parList,$spar);
#          push(@parList,$epar);
#    }
#
#    return @parList;
#}
#
#sub partition_inst_of_pin {
#    my $pin = shift;
#    #my ($up_pin_hier, $up_module, $sub_inst, $sub_ref, $sub_ref_pin) = get_pin_context_whier($pin);
#    my @up_pin_hier = get_pin_context_whier($pin);
#    $up_pin_hier = @up_pin_hier[0];
#    #@par_hier=grep(is_partition_module($_->[1]),get_hier_list_txt (-of_pin => $pin)); 
#    #@chiplet_hier=grep(is_chiplet_module($_->[1]),get_hier_list_txt (-of_pin => $pin)); 
#    #$up_pin_hier = lc($par_hier[0]->[1]);
#    #$up_chiplet_hier = lc($chiplet_hier[0]->[1]);
#    if (!is_partition_inst($up_pin_hier)) {
#        my $uphier = get_up_inst ($up_pin_hier);
#        $up_pin_hier = $uphier;
#    }
#    return ($up_pin_hier);
#}
#
#
#
# merged from get_path_feeds.mpl 
# may update some functions

#sub get_subchiplet_port_infor {
#    # extract sub-chiplet port information and output [sub-chiplet name,$port_x,$port_y]
#    my ($vio) = @_;
#    my @inter_pin_arr = attr_of_vio(inter_pin_array,$vio);
#    my @temp;
#    my @rtn_list;
#    for(my $j =0;$j<scalar(@inter_pin_arr);$j = $j + 1) {
#        my $inter_pin = $inter_pin_arr[$j];
#        my @hier_list = get_hier_list_txt (-quiet,-of_net => $inter_pin);
#        foreach my $hier_list (@hier_list) {
#            if (is_chiplet_module($hier_list->[1])) {
#                my $base_name = $inter_pin;
#                my $pattern = $hier_list->[0];
#                $pattern =~ s/\[/\\\[/g;
#                $pattern =~ s/\]/\\\]/g;
#                $base_name =~ s/$pattern//g;
#                $base_name =~ s/\/$//;
#                my @pin_xy = get_pin_xy($inter_pin);
#                push(@temp,[$base_name,$pin_xy[0],$pin_xy[1]]);
#            }
#        }
#    }
#    push(@rtn_list,[$temp[0]->[0],$temp[0]->[1],$temp[0]->[2]]);
#    for(my $i = 1;$i < scalar(@temp);$i = $i + 1) {
#        my $i_pre = $i -1;
#        if (($temp[$i]->[1] != $temp[$i_pre]->[1]) or ($temp[$i]->[2] != $temp[$i_pre]->[2])) {
#            push(@rtn_list,[$temp[$i]->[0],$temp[$i]->[1],$temp[$i]->[2]]);
#        }
#    }
#    return @rtn_list;
#}
#
sub get_port_partition {
    # detect ports entry/leaving partition by location
    my ($port) = @_;
    my @all_partitions = get_all_par_insts () ;
    #if ($#main::all_partitions == -1) {
    #    @main::all_partitions = get_all_par_insts () ;
    #}
    #my ($useless_var,@all_partitions) = expand_chiplet_in_partitions("",());
    my @port_xy = get_pin_xy($port);

    foreach my $s_par (@all_partitions) {
        if (is_point_in_par($port_xy[0],$port_xy[1],$s_par)) {
            return $s_par;
        }
    }
    return "";
}
#
#sub is_pin_in_sub_chiplet {
#    my ($pin,@all_partitions) = @_;
#    # if pin is port, recognise with the pin-xy
#    # if pin is not port , recognise by hier name
#    if (not is_port($pin)) {
#        if (grep(is_chiplet_module($_->[1]),get_hier_list_txt (-quiet,-of_pin => $pin))) {
#            return 1;
#        } else {
#            return 0;
#        }
#    } else {
#        my $par = get_port_partition($pin,@all_partitions);
#        if (attr_of_ref(is_chiplet,attr_of_cell(ref_name,$par))) {
#            return 1;
#        } else {
#            return 0;
#        }
#    }
#}
#
#
#
##----------------------------------------------------------#
## get_vio_feeds(path_of_id(408:35793));
## path_type: par2sub,sub2par,sub2sub,par2par
## start -> sub-chiplet port -> sub-chiplet end
## sub-chiplet start -> sub-chiplet port -> end 
## keep using tl_get_crossed_edges outside sub-chiplet
## use draw line to calculate feedthrough inside sub-chiplet
##----------------------------------------------------------#
#
#sub get_vio_feeds {
#    my ($vio,@all_partitions) = @_;
#    my @inter_pin_arr = attr_of_vio(inter_pin_array,$vio);
#    my $start_pin = name_of_pin(attr_of_vio(start_pin,$vio));
#    $start_pin =~ s/checkpin.*//g;
#    my $end_pin = name_of_pin(attr_of_vio(end_pin,$vio));
#
#    my @start_xy= get_pin_xy($start_pin);
#    my @end_xy= get_pin_xy($end_pin);
#    my @sub_chiplet_port_infor = ();
#    my $start_par;
#    my $end_par;
#    if (not is_port($start_pin)) {
#        $start_par = partition_inst_of_pin($start_pin);
#    } else {
#        $start_par = get_port_partition($start_pin,@all_partitions);
#    }
#    if (not is_port($end_pin)) {
#        $end_par = partition_inst_of_pin($end_pin);
#    } else {
#        $end_par = get_port_partition($end_pin,@all_partitions);
#    }
#
#    my @feed_pars = ();
#    if (is_pin_in_sub_chiplet($start_pin,@all_partitions) and is_pin_in_sub_chiplet($end_pin,@all_partitions)) {
#       # pathe start and end in sub-chiplet
#       # draw fly line from start point, through ports(if any), to end point
#       # calculate all the through partitions
#       push(@sub_chiplet_port_infor,[$start_par,$start_xy[0],$start_xy[1]]);
#       foreach my $port_infor (get_subchiplet_port_infor($vio)) {
#           push(@sub_chiplet_port_infor,$port_infor);
#       }
#       push(@sub_chiplet_port_infor,[$end_par,$end_xy[0],$end_xy[1]]);
#       my @feed_par_sub2sub;
#       for (my $i = 1; $i < scalar(@sub_chiplet_port_infor);$i = $i + 1) {
#           my @temp = calculate_feed_pars_by_line($sub_chiplet_port_infor[$i-1]->[1],$sub_chiplet_port_infor[$i-1]->[2],$sub_chiplet_port_infor[$i]->[1],$sub_chiplet_port_infor[$i]->[2],@all_partitions);
#           foreach my $s_item (@temp) {
#               push(@feed_par_sub2sub,$s_item);
#           }
#       }
#       @feed_pars = @feed_par_sub2sub;
#    } elsif (is_pin_in_sub_chiplet($start_pin,@all_partitions) and not is_pin_in_sub_chiplet($end_pin,@all_partitions)) {
#       # paths start in sub-chiplet and ends outside sub-chiplet
#       # A: use feedthrough infer by tl_get_crossed_edges from sub-chiplet port to endpoint 
#       # B: draw fly line from start to sub-chiplet port and calculate through partitions
#       # replace sub-chiplet name in A with the leaving partition in B
#       push(@sub_chiplet_port_infor,[$start_par,$start_xy[0],$start_xy[1]]);
#       foreach my $port_infor (get_subchiplet_port_infor($vio)) {
#           push(@sub_chiplet_port_infor,$port_infor);
#       }
#       my @feed_par_sub2par = tl_get_crossed_edges(-quiet,-min_cross => $sub_chiplet_port_infor[-1]->[1],$sub_chiplet_port_infor[-1]->[2],$end_xy[0],$end_xy[1],$sub_chiplet_port_infor[-1]->[0],$end_par);
#       my @feed_par_sub2sub;
#       for (my $i = 1 ; $i < scalar(@sub_chiplet_port_infor); $i = $i + 1) {
#           my @temp = calculate_feed_pars_by_line($sub_chiplet_port_infor[$i-1]->[1],$sub_chiplet_port_infor[$i-1]->[2],$sub_chiplet_port_infor[$i]->[1],$sub_chiplet_port_infor[$i]->[2],@all_partitions);
#           foreach my $s_item (@temp) {
#               push(@feed_par_sub2sub,$s_item);
#           }
#       }
#       foreach my $item (@feed_par_sub2sub) {
#           push(@feed_pars,$item);
#       }
#       for (my $i = 1;$i < scalar(@feed_par_sub2par);$i = $i + 1) {
#           push(@feed_pars,$feed_par_sub2par[$i]->[0]);
#       }
#    } elsif (not is_pin_in_sub_chiplet($start_pin,@all_partitions) and is_pin_in_sub_chiplet($end_pin,@all_partitions)) {
#       # paths start outside sub-chiplet and ends inside sub-chiplet
#       # A: use feedthrough infer by tl_get_crossed_edges from startpoint to sub-chiplet port ()
#       # B: draw fly line from sub-chiplet port to endpoints and calculate through partitions
#       # replace sub-chiplet name in A with the entry partition in B
#       @sub_chiplet_port_infor = get_subchiplet_port_infor($vio);
#       push(@sub_chiplet_port_infor,[$end_par,$end_xy[0],$end_xy[1]]);
#       my @feed_par_par2sub = tl_get_crossed_edges(-quiet,-min_cross => $start_xy[0],$start_xy[1],$sub_chiplet_port_infor[0]->[1],$sub_chiplet_port_infor[0]->[2],$start_par,$sub_chiplet_port_infor[0]->[0]);
#       my @feed_par_sub2sub;
#       for (my $i = 1; $i < scalar(@sub_chiplet_port_infor);$i = $i + 1) {
#           my @temp = calculate_feed_pars_by_line($sub_chiplet_port_infor[$i-1]->[1],$sub_chiplet_port_infor[$i-1]->[2],$sub_chiplet_port_infor[$i]->[1],$sub_chiplet_port_infor[$i]->[2],@all_partitions);
#           foreach my $s_item (@temp) {
#               push(@feed_par_sub2sub,$s_item);
#           }
#       }
#       for (my $i = 0;$i < (scalar(@feed_par_par2sub) -1);$i =  $i + 1) {
#           push(@feed_pars,$feed_par_par2sub[$i]->[0]);
#       }
#       foreach my $s_item (@feed_par_sub2sub) {
#           push(@feed_pars,$s_item);
#       }
#    } else {
#       # paths start outside sub-chiplet and ends  outside sub-chiplet 
#       # put this branch here in case we want to change the structure and method for feedpars analysis in retime_retour check.
#       if ($start_par ne $end_par) {
#           my @par2par_feed_pars = tl_get_crossed_edges(-quiet,-min_cross => $start_xy[0],$start_xy[1],$end_xy[0],$end_xy[1],$start_par,$end_par);
#           foreach my $s_item (@par2par_feed_pars) {
#               push(@feed_pars,$s_item->[0]);
#           }
#       } else {
#           @feed_pars = $end_par;
#       }
#    }
#    # we will see feedpars like : Par1 Par2 SubChiplet/ParA SubChiplet/ParA SubChiplet/ParB, require uniq before output feedpars
#    @feed_pars = uniq_feedpars_in_arr(@feed_pars);
#    return @feed_pars;
#}
#
#
#sub calculate_feed_pars_by_line {
#    # this method draws fly line from start_point to end_point, and calculate which partitions belongs for all extend points 
#    # and extend the point with step 50 , diagonise each point on the flyline if it's in the specific partitions
#    my ($start_x,$start_y,$end_x,$end_y,@all_partitions) = @_;
#    my @feed_pars = ();
#    #my ($useless_var,@all_partitions) = expand_chiplet_in_partitions("",());
#
#    my $start_par;
#    my $end_par;  
#    my $cal_step = 100;  
#    foreach my $s_par (@all_partitions) {
#        if (is_point_in_par(($start_x),($start_y),$s_par)) {
#            $start_par = $s_par;
#        }
#    } 
#    foreach my $s_par (@all_partitions) {
#        if (is_point_in_par(($end_x),($end_y),$s_par)) {
#            $end_par = $s_par;
#        }
#    } 
#
#    if ($s_par ne $end_par) {
#        push(@feed_pars,$start_par);
#        if ($start_x != $end_x) {
#            my $k = ($start_y -$end_y)/($start_x -$end_x);
#            my $b = ($start_x*$end_y - $start_y*$end_x)/($start_x -$end_x);
#            # detect startpoint and endpoint partition, it could be tricky if not move the corrdinate aside
#            # move corrdinate aside
#            if ($start_x < $end_x and $start_y < $end_y) {
#                $start_x = $start_x + 1;
#                $start_y = $start_y + 1;
#                $end_x = $end_x -1;
#                $end_y = $end_y -1;
#            } elsif ($start_x < $end_x and $start_y > $end_y) {
#                $start_x = $start_x + 1;
#                $start_y = $start_y - 1;
#                $end_x = $end_x - 1;
#                $end_y = $end_y + 1;
#            } elsif ($start_x > $end_x and $start_y > $end_y) {
#                $start_x = $start_x - 1;
#                $start_y = $start_y - 1;
#                $end_x = $end_x + 1;
#                $end_y = $end_y + 1;
#            } else {
#                $start_x = $start_x - 1;
#                $start_y = $start_y + 1;
#                $end_x = $end_x + 1;
#                $end_y = $end_y - 1;
#            }
#            if ($start_x < $end_x) {
#                for(my $point_x = $start_x + $cal_step;$point_x < $end_x;$point_x = $point_x + $cal_step) {
#                    my $point_x = int($point_x);
#                    my $point_y = int($point_x*$k + $b);
#                    my $hit_par_inside_chiplet = 0;
#                    for my $s_par (@all_partitions) {
#                        if (is_point_in_par($point_x,$point_y,$s_par)) {
#                            push(@feed_pars,$s_par);
#                            $hit_par_inside_chiplet = 1;
#                        }
#                    }
#                    if ($hit_par_inside_chiplet == 0) {
#                        push(@feed_pars,"NULL");
#                    }
#                }
#            } else {
#                foreach my $s_par (@all_partitions) {
#                    if (is_point_in_par(($start_x - 1),($start_y - 1),$s_par)) {
#                        push(@feed_pars,$s_par);
#                    }
#                }
#                for(my $point_x = $start_x - $cal_step;$point_x > $end_x;$point_x = $point_x - $cal_step) {
#                    my $point_x = int($point_x);
#                    my $point_y = int($point_x*$k + $b);
#                    my $hit_par_inside_chiplet = 0;
#                    for my $s_par (@all_partitions) {
#                        if (is_point_in_par($point_x,$point_y,$s_par)) {
#                            push(@feed_pars,$s_par);
#                            $hit_par_inside_chiplet = 1;
#                        }
#                    }
#                    if ($hit_par_inside_chiplet == 0) {
#                        push(@feed_pars,"NULL");
#                    }
#                }
#            }
#            push(@feed_pars,$end_par);
#        } else {
#           my $point_x = $start_x;
#           if ($start_y > $end_y) {
#               for(my $point_y = $start_y - $cal_step;$point_y > $end_y;$point_y = $point_y -$cal_step) {
#                    my $hit_par_inside_chiplet = 0;
#                    for my $s_par (@all_partitions) {
#                        if (is_point_in_par($point_x,$point_y,$s_par)) {
#                            push(@feed_pars,$s_par);
#                            $hit_par_inside_chiplet = 1;
#                        }
#                    }
#                    if ($hit_par_inside_chiplet == 0) {
#                        push(@feed_pars,"NULL");
#                    }
#               }
#           } else {
#               for(my $point_y = $start_y + $cal_step;$point_y < $end_y;$point_y = $point_y + $cal_step) {
#                    my $hit_par_inside_chiplet = 0;
#                    for my $s_par (@all_partitions) {
#                        if (is_point_in_par($point_x,$point_y,$s_par)) {
#                            push(@feed_pars,$s_par);
#                            $hit_par_inside_chiplet = 1;
#                        }
#                    }
#                    if ($hit_par_inside_chiplet == 0) {
#                        push(@feed_pars,"NULL");
#                    }
#               }
#           }
#        }
#        push(@feed_pars,$end_par);
#    } else {
#        push (@feed_pars,@start_par);
#    }
#    @feed_pars = uniq_feedpars_in_arr(@feed_pars);
#    return @feed_pars;
#}
#
#sub uniq_feedpars_in_arr {
#    # uniq A A B B B C C -> A B C
#    my (@feed_pars) = @_;
#    my @feed_pars_uniq = ($feed_pars[0]);
#    for (my $i = 1;$i<scalar(@feed_pars);$i = $i + 1) {
#        if ($feed_pars[$i -1] ne $feed_pars[$i]) {
#            push(@feed_pars_uniq,$feed_pars[$i]);
#        }
#    }
#    return @feed_pars_uniq;
#}


