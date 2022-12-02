
#' Read India Data Sets
#' 
#' @examples 
#' 
#' suppressMessages({
#'   dat = india_data("2020-04-01", "2020-06-01")
#'   dat = india_data("2020-04-01", "2020-07-01")
#' })
#' 
#' 
#' @importFrom readr read_csv
#' @importFrom magrittr `%>%`
#' @importFrom dplyr select filter arrange mutate tibble left_join n
india_data = function(date_initial, date_final) {
    # Set date ranges ----
    # date_initial = as.Date("2020-04-01")
    # date_final = as.Date("2020-06-01")

    # date_initial = as.Date("2020-04-01")
    # date_final = as.Date("2020-07-01")

    # Import data ----
    # confirmed
    case_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
    confirmed = read_csv(case_url)

    # deaths
    death_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
    deaths = read_csv(death_url)

    # recoveries
    recovery_url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv"
    recoveries = read_csv(recovery_url)

    # Process and merge data ----
    india_confirmed = confirmed %>%
    filter(`Country/Region` == "India") %>% #, is.na(`Province/State`)
    select(-`Province/State`, -Lat, -Long, -`Country/Region`)
    india_confirmed = tibble(date = colnames(india_confirmed),
                            cases_total = as.numeric(india_confirmed[1, ]))

    india_deaths = deaths %>%
    filter(`Country/Region` == "India") %>% #, is.na(`Province/State`)
    select(-`Province/State`, -Lat, -Long, -`Country/Region`)
    india_deaths = tibble(date = colnames(india_deaths),
                            deaths_total = as.numeric(india_deaths[1, ]))

    india_recovered = recoveries %>%
    filter(`Country/Region` == "India") %>% #, is.na(`Province/State`)
    select(-`Province/State`, -Lat, -Long, -`Country/Region`)
    india_recovered = tibble(date = colnames(india_recovered),
                            recoveries_total = as.numeric(india_recovered[1, ]))

    data = india_confirmed %>% left_join(india_recovered, by = "date")
    data = data %>% left_join(india_deaths, by = "date")

    data =
    data %>%
    mutate(date = as.Date(date, format = "%m/%d/%y"))

    data =
    data %>%
    mutate(total_removed = deaths_total + recoveries_total,
            active_cases = cases_total - total_removed,
            day = 1:n(),
            I = active_cases,
            R = total_removed)

    data =
    data %>%
    filter(I < 2e7, R > -3e5)

    data =
    data %>%
    filter(date >= date_initial, date <= date_final)

    data = 
    data %>%
    select(date, I, R)

    return(data)
}