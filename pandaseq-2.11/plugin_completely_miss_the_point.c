#include<errno.h>
#include<stdlib.h>
#include<pandaseq-plugin.h>

HELP("Filter out sequences that have mismatches in the overlap region.", "completely_miss_the_point:mismatches");

VER_INFO("1.0");

static bool check_func(
	PandaLogProxy logger,
	const panda_result_seq *sequence,
	void *user_data) {
	(void) logger;

	return sequence->overlap_mismatches <= (size_t) *(int *) user_data;
}

OPEN {
	int mismatches;

	(void) precheck;

	if (args == NULL || *args == '\0') {
		panda_log_proxy_write_str(logger, "Please supply the maximum allowed mismatches.\n");
		return false;
	}
	errno = 0;
	mismatches = strtol(args, NULL, 10);
	if (errno != 0 || mismatches < 0 || (size_t) mismatches > PANDA_MAX_LEN) {
		panda_log_proxy_write_str(logger, "Bad maximum allowed mismatches.\n");
		return false;
	}
	*check = check_func;
	*user_data = PANDA_STRUCT_DUP(&mismatches);
	*destroy = free;
	return true;
}
