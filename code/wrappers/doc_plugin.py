def doc_plugin(subroutine_line, name, run_type=None):

    table_string = []

    l = 0

    line = subroutine_line[l]

    while line.startswith("_COMMENT"):

        table_string.append(line.replace("_COMMENT>", ""))

        l += 1

        line = subroutine_line[l]

    return table_string
