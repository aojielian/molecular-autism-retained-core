# 00_session_info.R
# Export session info from HPC environment

sink("env/sessionInfo.txt")
sessionInfo()
sink()
