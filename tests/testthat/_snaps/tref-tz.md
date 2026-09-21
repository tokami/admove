# print.admove_tref shows origin, units and period

    Code
      print(create_tref(origin = as.Date("2025-01-01"), units = "month"))
    Output
      admove time reference
        origin:         2025-01-01 UTC
        units:          month
        period:         12
      

---

    Code
      print(create_tref())
    Output
      admove time reference
        origin:         not specified
        units:          not specified
        period:         none
      

