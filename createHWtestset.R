library(dplyr)

make_obs_hws <- function(city, print_city = TRUE){
  
  if(print_city) print(city)
  
  datafr <- readRDS(paste0("~/tmp/NMMAPS/", city, ".rds")) %>%
    dplyr::select(death, accident, date, agecat, tmpd, dptp, dow) %>%
    dplyr::mutate(death = death + accident)
  
  threshold <- quantile(datafr$tmpd[datafr$agecat == "under65"],
                        0.98, na.rm = TRUE)
  df_hws <- datafr %>%
    dplyr::filter(agecat == "under65") %>%
    select(date, tmpd) %>%
    futureheatwaves::IDHeatwavesCPPwrapper(threshold = threshold, numDays = 2)
  dist.tmpd <- ecdf(df_hws$tmpd)
  
  datafr <- datafr %>%
    dplyr::left_join(select(df_hws, -tmpd), by = "date")
  
  hw.list <- datafr %>%
    dplyr::filter(hw == 1) %>%
    dplyr::filter(agecat == "under65") %>%
    dplyr::group_by(hw.number) %>%
    dplyr::summarize(mean.temp = mean(tmpd),
                     max.temp = max(tmpd),
                     min.temp = min(tmpd),
                     length = length(unique(date)),
                     start.date = first(date),
                     end.date = last(date),
                     start.doy = yday(first(date)),
                     start.month = month(first(date)),
                     days.above.80 = sum(tmpd > 80),
                     days.above.85 = sum(tmpd > 85),
                     days.above.90 = sum(tmpd > 90),
                     days.above.95 = sum(tmpd > 95),
                     days.above.99th = sum(tmpd > quantile(datafr$tmpd, .99,
                                                           na.rm = TRUE)),
                     days.above.99.5th = sum(tmpd > quantile(datafr$tmpd, .995,
                                                             na.rm = TRUE))) %>%
    dplyr::mutate(first.in.season = ifelse(year(start.date) != 
                                             year(lag(start.date)),
                                           1, 0),
                  first.in.season = ifelse(is.na(first.in.season), 1, 
                                           first.in.season),
                  mean.temp.quantile = dist.tmpd(mean.temp),
                  max.temp.quantile = dist.tmpd(max.temp),
                  min.temp.quantile = dist.tmpd(min.temp),
                  mean.temp.1 = mean(datafr$tmpd, na.rm = TRUE),
                  mean.summer.temp = mean(datafr$tmpd[month(datafr$date) %in% 5:9], 
                                          na.rm = TRUE))
  
  datafr <- datafr %>%
    dplyr::mutate(time = scale(as.numeric(date), scale = FALSE, center = TRUE),
                  hw.number = factor(hw.number))
  city.mod <- glm(death ~ agecat + dow + splines::ns(time, 7 * 19) +
                    hw.number,
                  family = quasi(link = log, variance = mu),
                  data = datafr,
                  control=glm.control(epsilon=10E-8, maxit = 10000))
  
  hw.coef.rows <- grep("hw.number", names(coef(city.mod)))
  hw.coefs <- summary(city.mod)$coef[hw.coef.rows, 1:2] %>%
    dplyr::tbl_df() 
  hw.coefs <- dplyr::mutate(hw.coefs, hw.number = row.names(hw.coefs),
                  hw.number = gsub("hw.number", "", hw.number),
                  hw.number = as.numeric(hw.number))

  hw.list <- hw.list %>%
    ungroup() %>%
    left_join(hw.coefs, by = "hw.number") %>%
    mutate(city = city)
  
  return(hw.list)
}

city_list <- readr::read_csv("cities.csv")


hw_data <- lapply(city_list$city, make_obs_hws)
hw_data <- do.call("rbind", hw_data)

## Some checks
# hw_data %>% 
#   filter(max.temp.quantile < 0.995 & days.above.99.5th > 0) %>%
#   nrow()
# hw_data %>% 
#   filter(max.temp.quantile < 0.99 & days.above.99th > 0) %>%
#   nrow()
# hw_data %>% filter(mean.temp > max.temp) %>% nrow()
# hw_data %>% filter(min.temp > mean.temp) %>% nrow()
# hw_data %>% filter(days.above.99.5th > days.above.99th) %>% nrow()
# hw_data %>%
#   mutate(start.date = lubridate::ymd(start.date)) %>%
#   filter(start.month != lubridate::month(start.date)) %>%
#   nrow()

save(hw_data, file = "hwdata.RData")
load("hwdata.RData")

# Add population values and land area
pops <- readr::read_csv("nmmaps_counties.csv") %>%
  dplyr::select(city, pop) %>%
  left_join(readr::read_csv("land_area.csv"), by = "city") %>%
  group_by(city) %>%
  summarize(pop100 = sum(pop),
            arealand = sum(arealand),
            pop.density = pop100 / arealand)

hw_data <- hw_data %>% 
  left_join(pops, by = "city")

save(hw_data, file = "hwdatapop.RData")
load("hwdatapop.RData")

## For each heat wave, get posterior estimate and standard error
## after pooling all estimates
library(tlnise)
set.seed(21)
seed <- round(10000 * runif(1))

hw.pooled <- tlnise(hw_data$Estimate, hw_data$"Std. Error"^2,
                    seed = seed, maxiter = 5000, prnt = FALSE)

hw_data$post.estimate <- hw.pooled$theta
hw_data$post.se <- hw.pooled$SDtheta

save(hw_data, file = "data-raw/ListOfAllHeatwaves.Rdata")

# # Checks on effects of using posteriors
# library(ggplot2)
# y.range <- range(exp(hw_data$Estimate), exp(hw_data$post.estimate))
# library(scales)
# ggplot(hw_data, aes(x = pop100, y = exp(Estimate))) + 
#   geom_point(alpha = 0.3) + 
#   ylab("Relative risk") + xlab("Population") + 
#   scale_y_continuous(limits = y.range) + 
#   scale_x_continuous(labels = comma)
# ggplot(hw_data, aes(x = pop100, y = exp(post.estimate))) + 
#   geom_point(alpha = 0.3) + 
#   ylab("Relative risk") + xlab("Population") + 
#   scale_y_continuous(limits = y.range) + 
#   scale_x_continuous(labels = comma)
# 
# 


