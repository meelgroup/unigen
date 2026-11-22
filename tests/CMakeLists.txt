# unigen exits with an error if an extended sample is not a solution of the input
set(CNF ${CMAKE_CURRENT_SOURCE_DIR}/cnf)
foreach(f arjun-bve-def-1 arjun-bve-def-2 no-projection)
    add_test(NAME ${f} COMMAND unigen-bin --verb 0 --samples 200 --fullsol 1 ${CNF}/${f}.cnf)
    add_test(NAME ${f}-eappmc0 COMMAND unigen-bin --verb 0 --samples 200 --fullsol 1 --eappmc 0 ${CNF}/${f}.cnf)
    add_test(NAME ${f}-noarjun COMMAND unigen-bin --verb 0 --samples 200 --arjun 0 ${CNF}/${f}.cnf)
endforeach()
add_test(NAME unsat COMMAND unigen-bin --verb 0 ${CNF}/unsat.cnf)
set_tests_properties(unsat PROPERTIES WILL_FAIL TRUE)
