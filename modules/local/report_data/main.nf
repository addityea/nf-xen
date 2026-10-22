process REPORT_DATA {
    label 'process_single'

    container params.container_python

    input:
    path(figures)
    path(summaries)
    path params_file

    output:
    path("report_data.json"), emit: data
    path("report.md"),        emit: markdown
    path("figures/*"),        emit: figures

    script:
    template 'merge_report_data.py'
}
