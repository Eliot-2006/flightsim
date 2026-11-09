# Strict warnings + deterministic FP for our own targets only.
# Sanitizers are applied per-target (not globally) so third-party configure
# probes (e.g. Google Benchmark's regex check) are unaffected.
option(FLIGHTSIM_SANITIZE "Build our targets with ASan+UBSan" OFF)

function(flightsim_target_options target)
  target_compile_options(${target} PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Werror -ffp-contract=off)
  if(FLIGHTSIM_SANITIZE)
    target_compile_options(${target} PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(${target} PRIVATE -fsanitize=address,undefined)
  endif()
endfunction()
