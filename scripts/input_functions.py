import pdb

def analysis_df_value(wildcards, analysis_df, column_name):
	
	# get a value from the row of the df corresponding to this analysis condition
	unique = f"{wildcards.dset}"

	return analysis_df.loc[(analysis_df['analysis_condition'] == unique).idxmax(), column_name] 


def get_input_reads(wildcards, read_num):
	assert read_num in (1, 2)
	if read_num == 1:
		return f"{analysis_df_value(wildcards, 'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, 'R1_suffix')}"
	if read_num == 2:
		return f"{analysis_df_value(wildcards, 'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, 'R2_suffix')}"	


def get_polyidus_reads(wildcards, read_num):

	if analysis_df_value(wildcards, 'trim') == 1:
		if read_num == 1:
			return rules.trim.output.proc_r1
		else:
			return rules.trim.output.proc_r2
	else:
		if read_num == 1:
			return get_input_reads(wildcards, 1)
		else:
			return get_input_reads(wildcards, 2)
			
def get_vifi_resource(wildcards, analysis_df, resource_name):
	"""Get resources required for vifi"""
	host_idx = analysis_df[(analysis_df['host'] == wildcards.host) & (analysis_df['tool'] == 'vifi')].index[0]
	return analysis_df.loc[host_idx, resource_name]
	
def resources_list_with_min_and_max(file_name_list, attempt, mult_factor=2, minimum = 100, maximum = 50000):
	
	resource = int(sum([os.stat(file).st_size/1e6 for file in file_name_list]))
	
	resource = min(maximum, resource)
	
	return max(minimum, resource)

