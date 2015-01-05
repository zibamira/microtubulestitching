# Automated tests for hxalignmicrotubules pass

 - Purpose: Run automated tests for hxalignmicrotubules and dependencies.
 - Preconditions: none.
 - Platforms: all.

## Steps

Use gtest to run all tests of hxalignmicrotubules and dependencies.  But disable
slow hxcore tests:

    ./product/bin/zibamira -test \
        --packages=hxalignmicrotubules,pgmpointmatching \
        --gtest_filter=*-*_E7MS*:*HxWorkAreaTest*:*QxDevWizardTest*:*QxHistogramWorkerTest*

All tests should pass.
