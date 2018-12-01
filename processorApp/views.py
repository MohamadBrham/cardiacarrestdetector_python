from django.shortcuts import render

from django.http import HttpResponse
from . import PreProcessor


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")
# Create your views here.


def result(request):
    prediction, peaks= PreProcessor.classify(request.POST['data2'])
    # m = request.POST['data']
    return HttpResponse(prediction)
