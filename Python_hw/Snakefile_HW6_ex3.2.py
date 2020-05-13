def other_function ( * args , ** kwargs ) :
				return ""

def some_function ( * args , ** kwargs ) :
				return ""
@workflow.rule(name='some', lineno=6, snakefile='/home/aleksandar/snakemake_tutorial/examples/example_03/example_03.2.smk')

@workflow.output( "sorted_reads/{sample,[A-Za-z0-9]+}.bam"
)

@workflow.input(
				"plot.svg" ,
				some_function ( "condition1" ) ,
				[ "stats1.svg" , "stats2.svg" ] ,
				peaks1 = "cond1.bed" ,
				peaks2 = "cond2.bed" ,
				stats = other_function ( "condition1" )

)

@workflow.params(
#config["paths"]["{sample}"] # wrong, 'sample' isn't defined yet
				lambda wildcards : config [ "paths" ] [ wildcards . sample ] # correct
)
@workflow.norun()
@workflow.run
def __rule_some(input, output, params, wildcards, threads, resources, log, version, rule, conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, shadow_dir, edit_notebook):
	pass

