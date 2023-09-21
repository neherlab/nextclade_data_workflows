import datetime
import sys

def generate_dates(start_date_str, end_date_str):
    start_date = datetime.datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.datetime.strptime(end_date_str, "%Y-%m-%d")
    
    current_date = start_date
    while current_date <= end_date:
        print(f"designation_date\t{current_date.strftime('%Y-%m-%d')}")
        current_date += datetime.timedelta(days=1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python scriptname.py STARTDATE ENDDATE")
    else:
        generate_dates(sys.argv[1], sys.argv[2])
