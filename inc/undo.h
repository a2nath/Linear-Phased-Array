#pragma once

#include "spdlog/spdlog.h"

#define UNDOREDOSIZE 50

template<class T>
class UndoRedoBuffer {
private:
    struct Node {
        T settingData;
        Node* prev;
        Node* next;

        Node(T data) : settingData(data), prev(nullptr), next(nullptr) {}
    };
    Node* current; // Tracks the active state. Acts as your dynamic "head"
    Node* head;
    const size_t capacity;
    size_t size;

public:
    UndoRedoBuffer(size_t max_size = UNDOREDOSIZE) : current(nullptr), head(nullptr), capacity(max_size), size(0) {}

    // Destructor to clean up remaining memory
    ~UndoRedoBuffer() {

        // Delete forward
        while (head) {
            Node* nextNode = head->next;
            delete head;
            head = nextNode;
        }

    }

    inline void delete_first()
    {
        if (head)
        {
            Node* nextNode = head->next;
            delete head;
            head = nextNode;
        }
    }

    // 1. Send a request for a new setting
    void emplace(const T& newSetting) {
        if (size == capacity)
        {
            delete_first();
        }

        Node* newNode = new Node(newSetting);

        // Scenario A: First item in the buffer
        if (current == nullptr) {
            current = newNode;
            //printf("Added initial setting: " << current->settingData << "\n";
            return;
        }

        // Scenario B: User is in the middle of an undo chain.
        // We must delete all future "redo" nodes before adding the new one.
        Node* redoNode = current->next;
        while (redoNode != nullptr) {
            Node* temp = redoNode->next;
            //std::cout << "Deleting obsolete redo node: " << redoNode->settingData << "\n";
            delete redoNode;
            redoNode = temp;
            size--;
        }

        // Sever the link to the deleted future nodes
        current->next = nullptr;

        // Step C: Link the new node to the current node
        current->next = newNode;
        newNode->prev = current;

        // Move current pointer to the absolute latest setting
        current = newNode;
        //std::cout << "Added new setting: " << current->settingData << "\n";

        size++;
    }

    // 2. Undo operation
    void undo() {
        if (current && current->prev) {
            current = current->prev;
            //std::cout << "Undo successful. Current setting: " << current->settingData << "\n";
        }
        else {
            //std::cout << "Cannot Undo. You are at the oldest setting.\n";
        }
    }

    // 3. Redo operation
    void redo() {
        if (current && current->next) {
            current = current->next;
            //std::cout << "Redo successful. Current setting: " << current->settingData << "\n";
        }
        else {
            //std::cout << "Cannot Redo. You are at the latest setting.\n";
        }
    }

    bool has_changed() const
    {
        return current->prev->settingData == current->settingData;
    }

    // Helper to print the current active setting
    T& get_current() const {
        return current->settingData;
    }
};