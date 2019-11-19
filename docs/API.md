# Autometa User Registry Service (REST API)

## Usage

All responses will have the form

```json
{
    "data": "Mixed type holding the content of the response",
    "message": "Description of content"
}
```

Subsequent response definitions will only detail the expected value of the `data field`

### Register new user

**Definition**

`POST /user`

**Response**

- `200 OK` on success
- `409 Conflict` user with provided email already exists

```json
{
    "name": "Sam Waterworth",
    "email": "swaterworth@wisc.edu",
    "affiliation": "University of Wisconsin - Madison",
}
```

### List all projects

**Definition**

`GET /user/projects`

**Response**

- `200 OK` on success

```json
[
    {
        "identifier": "identifier",
        "name": "Stromatolite Study",
        "owner": "Sam Waterworth",
        "affiliation": "University of Wisconsin - Madison",
        "datasets": ["cape_recife", "schoenmakerscop"],
    },
    {
        "identifier": "identifier",
        "name": "Cranberry Study",
        "owner": "Kyle Wolf",
        "affiliation": "University of Wisconsin - Madison",
        "datasets": ["wild_cranberries", "wisconsin_bogs"],
    },
    {
        "identifier": "identifier",
        "name": "Fall 2020",
        "owner": "Tiny Earth Admin",
        "affiliation": "University of Wisconsin - Madison",
        "datasets": ["student1_soil_sample", "student2_soil_sample"],
    }
]
```

This would be response where `user = admin`, otherwise the response array would all share the same `owner` field.

### Registering a new project

**Definition**

`POST /user/projects`

**Arguments**

- `"name":string` a friendly name for this project
- `"owner":string` the owner/user of the project

**Response**

- `201 Created` on success
- `403 Forbidden` project with the given name already exists

```json
{
    "identifier": "identifier",
    "name": "Stromatolite Study",
    "owner": "Sam Waterworth",
    "affiliation": "University of Wisconsin - Madison",
    "datasets":[],
}
```

### Lookup project details

`GET /user/project/<identifier>`

**Response**

- `404 Not Found` if the project does not exist
- `200 OK` on success

```json
{
    "identifier": "identifier",
    "name": "Fall 2020",
    "owner": "Tiny Earth Admin",
    "affiliation": "University of Wisconsin - Madison",
    "datasets": ["picnic_point_soil", "picnic_point_firepit"],
}
```

### Delete a project

**Definition**

`DELETE /user/project/<identifier>`

**Response**

- `404 Not Found` if the project does not exist
- `204 No Content` on success
