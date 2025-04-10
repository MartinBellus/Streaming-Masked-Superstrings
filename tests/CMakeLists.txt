# Create symlink to data directory
set(DATA_DIR ${CMAKE_SOURCE_DIR}/data)
set(LINK_NAME ${CMAKE_BINARY_DIR}/data)

add_custom_command(
    OUTPUT ${LINK_NAME}
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${DATA_DIR} ${LINK_NAME}
    COMMENT "Creating symlink to data directory"
)

add_custom_target(create_data_symlink ALL DEPENDS ${LINK_NAME})

add_subdirectory(bloom_filter)
add_subdirectory(helper)
