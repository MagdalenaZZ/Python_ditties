


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct


def nosmall_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        #val = int(round(pct*total/100.0)+0.5)

        #print(pct,val)
        if pct < 5:
            val=' '
            pct=' '
            return (pct)

        return '{p:.1f}'.format(p=pct)
    return my_autopct



