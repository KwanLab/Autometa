from django.shortcuts import render

jobs = [
    {
        'user': 'schanana',
        'title': 'Job 1',
        'parameters': [1, 2, 3, 'four'],
        'date_run': 'August 27, 2018'
    },
    {
        'user': 'test',
        'title': 'Job 2',
        'parameters': [5, 6, 7, 'eight'],
        'date_run': 'November 29, 2018'
    }
]


def startpage(request):
    context = {'jobs': jobs}
    return render(request, 'startpage/home.html', context)


def about(request):
    return render(request, 'startpage/about.html', context={'title': 'about autometa'})
