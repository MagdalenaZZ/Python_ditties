
'''
#First, set a PROJECT_ID environment variable:
export PROJECT_ID=$(gcloud config get-value core/project)

#Next, create a new service account to access the Text-to-Speech API by using:
gcloud iam service-accounts create my-tts-sa  --display-name "my tts service account"

#Next, create credentials that your Python code will use to login as your new service account. Create and save these credentials as a ~/key.json JSON file by using the following command:
gcloud iam service-accounts keys create ~/key.json  --iam-account my-tts-sa@${PROJECT_ID}.iam.gserviceaccount.com

#Finally, set the GOOGLE_APPLICATION_CREDENTIALS environment variable, which is used by the Speech-to-Text client library, covered in the next step, to find your credentials. The environment variable should be set to the full path of the credentials JSON file you created:
export GOOGLE_APPLICATION_CREDENTIALS=~/key.json


webmd-324920

export GOOGLE_APPLICATION_CREDENTIALS='/Users/mz3/GAPI/webmd-324920-a268df3df41d.json'


'''

import google.cloud.texttospeech as tts

from google.oauth2 import service_account

credentials = service_account.Credentials.from_service_account_file('/Users/mz3/GAPI/webmd-324920-a268df3df41d.json')


def unique_languages_from_voices(voices):
    language_set = set()
    for voice in voices:
        for language_code in voice.language_codes:
            language_set.add(language_code)
    return language_set

# List languages
def list_languages():
    client = tts.TextToSpeechClient()
    response = client.list_voices()
    languages = unique_languages_from_voices(response.voices)
    print(f" Languages: {len(languages)} ".center(60, "-"))
    for i, language in enumerate(sorted(languages)):
        print(f"{language:>10}", end="\n" if i % 5 == 4 else "")

list_languages()

# List voices
def list_voices(language_code=None):
    client = tts.TextToSpeechClient()
    response = client.list_voices(language_code=language_code)
    voices = sorted(response.voices, key=lambda voice: voice.name)
    print(f" Voices: {len(voices)} ".center(60, "-"))
    for voice in voices:
        languages = ", ".join(voice.language_codes)
        name = voice.name
        gender = tts.SsmlVoiceGender(voice.ssml_gender).name
        rate = voice.natural_sample_rate_hertz
        print(f"{languages:<8} | {name:<24} | {gender:<8} | {rate:,} Hz")


def text_to_wav(voice_name: str, text: str):
    language_code = "-".join(voice_name.split("-")[:2])
    text_input = tts.SynthesisInput(text=text)
    voice_params = tts.VoiceSelectionParams(
        language_code=language_code, name=voice_name
    )
    audio_config = tts.AudioConfig(audio_encoding=tts.AudioEncoding.LINEAR16)
    client = tts.TextToSpeechClient()
    response = client.synthesize_speech(
        input=text_input, voice=voice_params, audio_config=audio_config
    )
    filename = f"{language_code}.wav"
    with open(filename, "wb") as out:
        out.write(response.audio_content)
        print(f'Generated speech saved to "{filename}"')

text_to_wav("en-AU-Wavenet-A", "Ulcerative colitis (UC) is type of inflammatory bowel disease that causes sores in the colon.")
text_to_wav("en-GB-Wavenet-B", "Ulcerative colitis (UC) is type of inflammatory bowel disease that causes sores in the colon.")
text_to_wav("en-IN-Wavenet-C", "Ulcerative colitis (UC) is type of inflammatory bowel disease that causes sores in the colon.")
text_to_wav("en-US-Wavenet-F", "Ulcerative colitis (UC) is type of inflammatory bowel disease that causes sores in the colon.")


def synthesize_text(text, title):
    """Synthesizes speech from the input string of text."""
    output=title + '.mp3'
    from google.cloud import texttospeech
    client = texttospeech.TextToSpeechClient()
    input_text = texttospeech.SynthesisInput(text=text)
    # Note: the voice can also be specified by name.
    # Names of voices can be retrieved with client.list_voices().
    voice = texttospeech.VoiceSelectionParams(
        language_code="en-GB",
        name="en-GB-Wavenet-B",
        ssml_gender=texttospeech.SsmlVoiceGender.MALE,
    )
    audio_config = texttospeech.AudioConfig(
        audio_encoding=texttospeech.AudioEncoding.MP3
    )
    response = client.synthesize_speech(
        request={"input": input_text, "voice": voice, "audio_config": audio_config}
    )
    # The response's audio_content is binary.
    with open(output, "wb") as out:
        out.write(response.audio_content)
        print('Audio content written to file ', output)

synthesize_text( "Ulcerative colitis (UC) is type of inflammatory bowel disease that causes sores in the colon.", "Ulc") 


