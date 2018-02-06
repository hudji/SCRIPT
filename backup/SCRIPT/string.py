s = "123123STRINGabcabc"
def find_between_r( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""

print find_between_r( s, "123", "abc" )