# GTA Housing Outlook | ECON 354
#
# Install once in RStudio:
#   install.packages(c("shiny", "plotly"))
#
# Put this script and the two CSVs in the same project folder. Set that
# folder as the working directory in RStudio, then run:
#   source("GTA_Housing_Outlook.R")
#
# You can also run Rscript GTA_Housing_Outlook.R from that folder.
# Edit these paths if your downloaded CSVs have a different suffix.

# The CSVs are read from the working directory if they are there. Otherwise the
# script falls back to the project folder below.
folder <- "/Users/mialevcovich/Desktop/fall courses 2026/ECON 354/econ354 term project"
if (file.exists("TRREB_GTA_Monthly_Market_Watch.csv")) folder <- getwd()
trreb_file <- file.path(folder, "TRREB_GTA_Monthly_Market_Watch.csv")
rate_file <- file.path(folder, "RATES_MORTGAGES-2013-01-01.csv")

test_start <- as.Date("2020-01-01")
test_end <- as.Date("2025-12-01")
scenario_origin <- as.Date("2026-08-01")
chart_origins <- seq(test_start, scenario_origin, by = "month")
model_labels <- c("Price-only forecast", "With listings")

# 1. Load the same inputs and select municipalities using pre-test sales.
market <- read.csv(trreb_file)
market <- subset(market, area_type == "municipality",
                 select = c(area, date, sales, new_listings, median_price))
market$date <- as.Date(market$date)

rate_header <- match('"OBSERVATIONS"', trimws(readLines(rate_file, warn = FALSE)))
if (is.na(rate_header)) stop("Could not find the Bank of Canada OBSERVATIONS table.")
rates <- read.csv(rate_file, skip = rate_header)
rates$date <- as.Date(rates$date)
stopifnot(!anyNA(market$date), !anyNA(rates$date),
          !anyDuplicated(market[c("area", "date")]), !anyDuplicated(rates$date))

pretest_sales <- aggregate(sales ~ area, subset(market, date < test_start), median)
municipalities <- sort(pretest_sales$area[pretest_sales$sales >= 100])
stopifnot(length(municipalities) == 15L)

# 2. Fit two direct six-month models separately for each municipality.
# The price-only model uses recent six-month growth and calendar month.
# The second model adds the ratio of total sales to total new listings
# over t-4 through t-2. Both models train on the same completed windows.
lag_months <- function(x, k) c(rep(NA, k), head(x, -k))
lead_months <- function(x, k) c(x[-seq_len(k)], rep(NA, k))
percent_growth <- function(log_growth) 100 * expm1(log_growth)

forecast_municipality <- function(local) {
  local <- local[order(local$date), ]
  stopifnot(all(local$date == seq(min(local$date), max(local$date), by = "month")),
            all(local$median_price > 0), all(local$sales >= 0),
            all(local$new_listings > 0))

  log_price <- log(local$median_price)
  local$recent_growth <- log_price - lag_months(log_price, 6)
  local$future_growth <- lead_months(log_price, 6) - log_price
  local$calendar_month <- factor(format(local$date, "%m"), levels = sprintf("%02d", 1:12))
  sales <- lag_months(local$sales, 2) + lag_months(local$sales, 3) +
    lag_months(local$sales, 4)
  listings <- lag_months(local$new_listings, 2) + lag_months(local$new_listings, 3) +
    lag_months(local$new_listings, 4)
  local$snlr <- sales / listings

  # A forecast made at t can only train on outcomes that end before t.
  target_dates <- seq(local$date[1], by = "month", length.out = nrow(local) + 6L)
  local$target_date <- target_dates[seq_len(nrow(local)) + 6L]
  complete <- complete.cases(local[c("recent_growth", "future_growth", "snlr")])
  origin_rows <- match(chart_origins, local$date)
  if (anyNA(origin_rows)) stop("Missing a forecast starting month for ", local$area[1])

  forecasts <- lapply(origin_rows, function(i) {
    current <- local[i, ]
    training <- local[complete & local$target_date < current$date, ]
    price_only <- lm(future_growth ~ recent_growth + calendar_month, data = training)
    with_listings <- lm(future_growth ~ recent_growth + calendar_month + snlr, data = training)
    data.frame(
      area = current$area,
      date = current$date,
      target_date = current$target_date,
      starting_price = current$median_price,
      actual = percent_growth(current$future_growth),
      price_only = percent_growth(unname(predict(price_only, current))),
      with_listings = percent_growth(unname(predict(with_listings, current)))
    )
  })
  do.call(rbind, forecasts)
}

market_by_area <- split(market[market$area %in% municipalities, ],
                        market$area[market$area %in% municipalities])
forecasts <- do.call(rbind, lapply(market_by_area, forecast_municipality))
forecasts$start_rate <- rates$V122667786[match(forecasts$date, rates$date)]
forecasts$end_rate <- rates$V122667786[match(forecasts$target_date, rates$date)]

# 3. Score the 72 original test starts. The 2026 extension is for display.
backtest <- subset(forecasts, date <= test_end)
stopifnot(nrow(backtest) == length(municipalities) * 72L,
          all(complete.cases(backtest[c("actual", "price_only", "with_listings")])),
          all(complete.cases(forecasts[c("price_only", "with_listings")])))

score <- function(prediction, actual) {
  error <- prediction - actual
  c(RMSE = sqrt(mean(error^2)), MAE = mean(abs(error)),
    Direction = 100 * mean(sign(prediction) == sign(actual)))
}

price_scores <- function(rows) {
  rbind(score(rows$price_only, rows$actual), score(rows$with_listings, rows$actual))
}

pooled_prices <- price_scores(backtest)
price_by_area <- lapply(split(backtest, backtest$area), price_scores)

# 4. Translate prices into mortgage payments using the same financing
# assumptions. For the backtest, both models get the same observed rates.
payment_factor <- function(rate_pct, years = 25) {
  monthly <- expm1(log1p(rate_pct / 200) / 6)
  result <- monthly / -expm1(-12 * years * log1p(monthly))
  result[which(monthly == 0)] <- 1 / (12 * years)
  result
}

payment_change <- function(price_growth, financing_ratio) {
  100 * ((1 + price_growth / 100) * financing_ratio - 1)
}

payments <- subset(backtest, is.finite(start_rate) & is.finite(end_rate))
payments$regime <- sign(payments$end_rate - payments$start_rate)
payment_scores <- function(rows) {
  financing <- payment_factor(rows$end_rate) / payment_factor(rows$start_rate)
  actual <- payment_change(rows$actual, financing)
  rbind(score(payment_change(rows$price_only, financing), actual),
        score(payment_change(rows$with_listings, financing), actual))
}

scenario <- function(row, future_rate, income_growth) {
  financing <- payment_factor(future_rate) / payment_factor(row$start_rate)
  income_ratio <- 1 + income_growth / 100
  starting_payment <- row$starting_price * 0.80 * payment_factor(row$start_rate)
  change <- payment_change(row$with_listings, financing)
  list(price = row$with_listings, payment = change,
       pti = 100 * ((1 + change / 100) / income_ratio - 1),
       break_even = 100 * (income_ratio / financing - 1),
       starting_payment = starting_payment,
       future_payment = starting_payment * (1 + change / 100))
}

# 5. Optional numerical check. This runs with base R only.
check_only <- identical(Sys.getenv("GTA_DASHBOARD_CHECK_ONLY"), "1")
if (check_only) {
  print(pooled_prices)
  print(payment_scores(payments))
  print(payment_scores(subset(payments, regime < 0)))
  print(price_by_area[["Oakville"]])
  print(subset(forecasts, area == "Oakville" & date >= as.Date("2026-04-01"),
               select = c(date, target_date, actual, price_only, with_listings)))
}

# 6. Render the interactive Shiny dashboard.
build_dashboard <- function() {
  library(shiny)
  library(plotly)

  ink <- "#10382e"
  green <- "#157458"
  rust <- "#ac5e38"
  blue <- "#627c8b"
  signed <- function(x) sprintf("%+.2f%%", x)
  dollars <- function(x) paste0("$", formatC(x, format = "f", digits = 0, big.mark = ","))
  month_label <- function(x) format(x, "%b %Y")
  arrow <- function(x) sprintf("%.2f to %.2f", x[1], x[2])
  wins <- sum(vapply(price_by_area, function(x) x[2, "RMSE"] < x[1, "RMSE"], logical(1)))
  falling <- payment_scores(subset(payments, regime < 0))
  plot_options <- function(chart) config(chart, displayModeBar = FALSE, responsive = TRUE)
  card <- function(...) tags$div(class = "card", ...)
  stat <- function(label, value, detail) {
    tags$div(class = "stat", tags$div(label), tags$b(value), tags$small(detail))
  }

  ui <- fluidPage(
    tags$head(tags$style(HTML("
      body {background:#f1f5f1;color:#10382e;font:16px/1.5 system-ui,sans-serif}
      main {max-width:1240px;margin:auto;padding:20px clamp(12px,4vw,36px) 42px}
      header {background:#10382e;color:#f7fbf4;padding:32px clamp(12px,4vw,36px);border-radius:12px}
      h1,h2 {font-family:Georgia,serif} h2 {margin:30px 0 8px} h3 {font-size:20px}
      .card,.stat {background:white;border:1px solid #d5e3d9;border-radius:12px;padding:18px;min-width:0}
      .stats,.scenario {display:grid;gap:14px;margin-top:14px}
      .stats {grid-template-columns:repeat(3,1fr)}
      .scenario {grid-template-columns:minmax(240px,1fr) minmax(0,1.35fr)}
      .stat b {display:block;font:28px Georgia,serif;margin:5px 0}
      .boundary {font-weight:650;color:#285e50}
      .boundary:before {content:'';display:inline-block;width:22px;height:12px;background:#e5f2eb;
        border-left:2px dashed #285e50;margin-right:8px}
      .detail {background:#eaf4ec;border-radius:7px;padding:10px;margin-top:10px}
      .control-row {display:flex;flex-wrap:wrap;gap:14px;justify-content:space-between}
      .table-scroll {overflow-x:auto} details {margin-top:18px}
      @media(max-width:700px) {.stats,.scenario {grid-template-columns:1fr}}
    "))),
    tags$main(
      tags$header(tags$div("ECON 354 · Greater Toronto Area"),
        h1("Resale prices and mortgage affordability"),
        p("Do sales and listings improve six-month price forecasts?")),
      tags$div(class = "stats",
        stat("Municipalities with lower price RMSE", paste(wins, "of", length(municipalities)),
             "After adding the listings ratio"),
        stat("GTA-wide price RMSE", arrow(pooled_prices[, "RMSE"]),
             "Percentage points of six-month growth"),
        stat("Payment RMSE when rates fell", arrow(falling[, "RMSE"]),
             "Using the same observed rates for both models")),
      tags$div(id = "forecast-results", class = "control-row",
        tags$div(h2("Where did listings help?"),
          p(paste("Forecast starting months through", month_label(scenario_origin),
                  "· Accuracy scores use the January 2020 to December 2025 backtest."))),
        selectInput("area", "Municipality", municipalities, "City of Toronto")),
      card(h3(textOutput("chart_area", inline = TRUE)),
        p("Six-month median resale price growth by forecast starting month."),
        uiOutput("boundary_key"),
        plotlyOutput("history", height = "400px"),
        p(class = "detail", textOutput("latest_forecast", inline = TRUE)),
        p(textOutput("boundary_note", inline = TRUE)),
        p(class = "detail", textOutput("area_accuracy", inline = TRUE)),
        p("RMSE gives more weight to large misses. MAE measures the average absolute miss.")),
      tags$div(style = "margin-top:14px", card(
        h3("Change in forecast error"), p("Price-only error minus error with listings."),
        radioButtons("metric", "Error measure", c("RMSE", "MAE"), inline = TRUE),
        plotlyOutput("rankings", height = "490px"))),
      h2("Mortgage affordability"),
      p(textOutput("scenario_window", inline = TRUE)),
      tags$div(class = "scenario",
        card(h3("Set your assumptions"), p(textOutput("starting_values", inline = TRUE)),
          sliderInput("future_rate", "Mortgage rate in six months (%)", 2, 8, 4.5, step = .05),
          sliderInput("income", "Household income growth over six months (%)", -8, 12, 2, step = .25),
          actionButton("reset", "Reset assumptions"),
          p("Illustrative purchase with 20% down and a 25-year amortization. Payments cover principal and interest.")),
        card(h3("Projected change over six months"), tableOutput("scenario_values"),
          plotlyOutput("break_even", height = "185px"),
          p(class = "detail", textOutput("scenario_note", inline = TRUE)),
          p(textOutput("payments_dollars", inline = TRUE)))),
      h2("Mortgage-payment backtest"),
      p("Both forecasts use the same observed rates and are compared with the payment implied by actual prices."),
      card(tags$div(class = "control-row",
          radioButtons("regime", "Mortgage rate environment",
            c("All months" = "all", "Rates fell" = "fall", "Rates rose" = "rise"), inline = TRUE),
          selectInput("scope", "Sample", c("All 15 municipalities" = "gta", "Selected municipality" = "local"))),
        p(textOutput("sample_note", inline = TRUE)),
        tags$div(class = "table-scroll", tableOutput("payment_accuracy"))),
      tags$details(class = "card", tags$summary("Data and methodology"),
        p("TRREB Market Watch provides monthly municipal median resale prices, sales and new listings. Changes in the mix of homes sold can affect the median."),
        p("Both direct models use recent six-month price growth and calendar-month indicators. The listings model adds the sales-to-new-listings ratio from months two through four before the forecast starts."),
        p("Mortgage rates use Bank of Canada series V122667786. The affordability scenario shows the relative change in the payment-to-income ratio under the chosen rate and income assumptions.")
      )
    )
  )

  server <- function(input, output, session) {
    current <- reactive(subset(forecasts, area == input$area))
    scenario_row <- reactive(subset(current(), date == scenario_origin))
    scenario_values <- reactive(scenario(scenario_row(), input$future_rate, input$income))
    output$chart_area <- renderText(input$area)

    forecast_boundary <- reactive({
      rows <- current()
      last_actual <- max(which(is.finite(rows$actual)))
      if (last_actual == nrow(rows)) return(NULL)
      rows[last_actual + 1L, ]
    })
    output$boundary_key <- renderUI({
      first <- forecast_boundary()
      if (!is.null(first)) {
        p(class = "boundary", paste("Forecast-only windows start in", month_label(first$date)))
      }
    })
    output$boundary_note <- renderText({
      first <- forecast_boundary()
      if (is.null(first)) return("")
      past <- subset(current(), is.finite(actual))
      paste("The first forecast-only window runs from", month_label(first$date), "to",
            paste0(month_label(first$target_date), "."), "Observed prices go through",
            paste0(month_label(max(past$target_date)), "."))
    })
    output$latest_forecast <- renderText({
      row <- scenario_row()
      paste(month_label(row$date), "to", month_label(row$target_date),
            "· With listings", signed(row$with_listings))
    })
    output$area_accuracy <- renderText({
      rmse <- price_by_area[[input$area]][, "RMSE"]
      sprintf("Backtest RMSE · Price-only %.2f pp · With listings %.2f pp", rmse[1], rmse[2])
    })

    output$history <- renderPlotly({
      rows <- current()
      window <- paste(month_label(rows$date), "to", month_label(rows$target_date))
      chart <- plot_ly() |>
        add_lines(x = rows$date, y = rows$actual, name = "Actual",
          line = list(color = blue, width = 2), text = window,
          hovertemplate = "%{text}<br>Actual %{y:+.2f}%<extra></extra>") |>
        add_lines(x = rows$date, y = rows$price_only, name = model_labels[1],
          line = list(color = rust, dash = "dash", width = 2), text = window,
          hovertemplate = "%{text}<br>Price-only %{y:+.2f}%<extra></extra>") |>
        add_lines(x = rows$date, y = rows$with_listings, name = model_labels[2],
          line = list(color = green, width = 3), text = window,
          hovertemplate = "%{text}<br>With listings %{y:+.2f}%<extra></extra>")
      first <- forecast_boundary()
      shapes <- list()
      if (!is.null(first)) {
        boundary <- as.character(first$date)
        shapes <- list(
          list(type = "rect", xref = "x", yref = "paper", x0 = boundary,
               x1 = as.character(max(rows$date)), y0 = 0, y1 = 1, layer = "below",
               fillcolor = "#e5f2eb", line = list(width = 0)),
          list(type = "line", xref = "x", yref = "paper", x0 = boundary, x1 = boundary,
               y0 = 0, y1 = 1, line = list(color = "#285e50", width = 2, dash = "dash"))
        )
      }
      plot_options(layout(chart, shapes = shapes, hovermode = "x",
        legend = list(orientation = "h", y = -.23),
        margin = list(l = 52, r = 24, t = 14, b = 86),
        xaxis = list(type = "date", title = "Forecast starting month", showgrid = FALSE,
                     range = as.character(range(rows$date)), tickformat = "%b %Y", nticks = 7),
        yaxis = list(title = "", ticksuffix = "%", gridcolor = "#e7eee9")))
    })

    output$rankings <- renderPlotly({
      gains <- vapply(price_by_area, function(scores) {
        scores[1, input$metric] - scores[2, input$metric]
      }, numeric(1))
      gains <- sort(gains)
      chart <- plot_ly(
        x = unname(gains), y = names(gains), key = names(gains), source = "rankings",
        type = "bar", orientation = "h",
        marker = list(color = ifelse(gains < 0, rust, green)),
        hovertemplate = "%{y}<br>%{x:+.2f} pp<extra></extra>"
      ) |>
        event_register("plotly_click") |>
        layout(xaxis = list(title = "Error reduction (percentage points)"),
          yaxis = list(title = "", categoryorder = "array", categoryarray = names(gains)),
          margin = list(l = 130, r = 24, t = 10, b = 56))
      plot_options(chart)
    })
    observeEvent(event_data("plotly_click", source = "rankings"), {
      clicked <- event_data("plotly_click", source = "rankings")
      if (length(clicked$key) == 1L && clicked$key %in% municipalities) {
        updateSelectInput(session, "area", selected = clicked$key)
      }
    })

    observeEvent(input$reset, {
      updateSliderInput(session, "future_rate", value = 4.5)
      updateSliderInput(session, "income", value = 2)
    })
    output$scenario_window <- renderText({
      row <- scenario_row()
      paste(row$area, "·", month_label(row$date), "to", month_label(row$target_date), "· Listings model")
    })
    output$starting_values <- renderText({
      row <- scenario_row()
      sprintf("Starting mortgage rate %.2f%% · Median price %s", row$start_rate, dollars(row$starting_price))
    })
    output$scenario_values <- renderTable({
      values <- scenario_values()
      data.frame(Measure = c("Median resale price", "Monthly payment", "Payment as a share of income"),
                 Change = signed(c(values$price, values$payment, values$pti)))
    }, striped = TRUE)
    output$break_even <- renderPlotly({
      values <- scenario_values()
      plot_ly(x = c(values$break_even, values$price),
        y = c("Break-even price growth", "Price forecast"), type = "bar", orientation = "h",
        marker = list(color = c("#8da8a9", green)),
        hovertemplate = "%{y}<br>%{x:+.2f}%<extra></extra>") |>
        layout(margin = list(l = 150, t = 10, r = 20, b = 34),
          xaxis = list(title = "", ticksuffix = "%"), yaxis = list(title = "")) |>
        plot_options()
    })
    output$scenario_note <- renderText({
      values <- scenario_values()
      burden <- if (abs(values$pti) < .005) "stays roughly unchanged" else if (values$pti < 0) "falls" else "rises"
      paste("With break-even price growth of", paste0(signed(values$break_even), ","),
            "the mortgage payment share of income", paste0(burden, "."))
    })
    output$payments_dollars <- renderText({
      values <- scenario_values()
      paste("Estimated monthly payment goes from", dollars(values$starting_payment),
            "to", paste0(dollars(values$future_payment), "."))
    })

    payment_sample <- reactive({
      rows <- payments
      if (input$scope == "local") rows <- subset(rows, area == input$area)
      if (input$regime == "fall") rows <- subset(rows, regime < 0)
      if (input$regime == "rise") rows <- subset(rows, regime > 0)
      rows
    })
    output$sample_note <- renderText({
      rows <- payment_sample()
      paste(length(unique(rows$date)), "forecast starting months ·",
            nrow(rows), "municipality-month observations")
    })
    output$payment_accuracy <- renderTable({
      rows <- payment_sample()
      validate(need(nrow(rows) > 0, "No paired mortgage rates for this sample."))
      scores <- payment_scores(rows)
      data.frame("Price forecast" = model_labels,
                 "Payment RMSE" = sprintf("%.2f pp", scores[, "RMSE"]),
                 "Payment MAE" = sprintf("%.2f pp", scores[, "MAE"]),
                 "Correct rise or fall" = sprintf("%.2f%%", scores[, "Direction"]),
                 check.names = FALSE)
    }, striped = TRUE)
  }
  shinyApp(ui, server)
}

if (!check_only) {
  packages <- c("shiny", "plotly")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install the dashboard packages first with install.packages(c('shiny', 'plotly')).",
         call. = FALSE)
  }
  shiny::runApp(build_dashboard(), launch.browser = TRUE)
}
