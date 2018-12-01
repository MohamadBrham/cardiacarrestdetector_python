from django.shortcuts import render

from django.http import HttpResponse, JsonResponse
from . import PreProcessor


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")
# Create your views here.




def result(request):
    data = request.POST['data']
    prediction, peaks= PreProcessor.classify(data)
    print(str(data))
    # m = request.POST['data']
    return JsonResponse({'Prediction': prediction})
