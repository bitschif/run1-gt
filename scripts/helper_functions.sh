#!/bin/bash
#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_warn() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [WARN] $1" >&2
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [ERROR] $1" >&2
}

check_file() {
    if [[ ! -f "$1" ]]; then
        log_error "File not found: $1"
        return 1
    fi
    return 0
}

check_tool() {
    if ! command -v "$1" &> /dev/null; then
        log_error "Tool not found: $1"
        return 1
    fi
    return 0
}

ensure_dir() {
    [[ !  -d "$1" ]] && mkdir -p "$1"
}

start_timer() {
    export STEP_START_TIME=$(date +%s)
}

end_timer() {
    local step_name="$1"
    local end_time=$(date +%s)
    local duration=$((end_time - STEP_START_TIME))
    log_info "${step_name} completed in ${duration} seconds"
    echo "${step_name},${duration}" >> "${LOG_DIR}/runtime.csv"
}

check_exit() {
    if [[ $?  -ne 0 ]]; then
        log_error "$1 failed"
        exit 1
    fi
    log_info "$1 completed"
}
