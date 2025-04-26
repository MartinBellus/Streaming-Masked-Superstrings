add_executable(fasta_test fasta_test.cpp)
add_dependencies(fasta_test create_data_symlink)
target_link_libraries(fasta_test io)
